# -*- coding: utf-8 -*-
# app.py — GHGSat C10 — 1 semana (v9e, duas colunas apenas)
# - Só GHGSAT-C10
# - Sem coluna de parâmetros (nenhuma UI no meio)
# - Lê TLE local ./tle/ghgsat.tle
# - Ao confirmar AOI, calcula semana (ASC/DESC) e plota footprints 5×5 km

import os, re, json, math
from datetime import datetime, timezone, timedelta

import streamlit as st
import folium
from streamlit_folium import st_folium
from shapely.geometry import Polygon, shape, mapping, box
from shapely.ops import unary_union
from shapely.affinity import rotate, translate
from pyproj import CRS, Transformer
from skyfield.api import EarthSatellite, load, wgs84
import numpy as np

st.set_page_config(page_title='GHGSat C10 — Semana (TLE offline)', layout='wide')

TLE_PATH = os.path.join(os.path.dirname(__file__), 'tle', 'ghgsat.tle')

# ---------- helpers geom ----------
def utm_crs_from_lonlat(lon: float, lat: float):
    from pyproj import CRS
    zone = int((lon + 180) // 6) + 1
    south = lat < 0
    return CRS.from_dict({'proj': 'utm', 'zone': zone, 'south': south})

def to_proj(geom_geojson, crs_to):
    from shapely.geometry import shape, Polygon
    from shapely.ops import unary_union
    from pyproj import Transformer
    geom = shape(geom_geojson)
    tr = Transformer.from_crs('epsg:4326', crs_to, always_xy=True)
    def P(pt): return tr.transform(pt[0], pt[1])
    if geom.geom_type == 'Polygon':
        ext = [P(c) for c in geom.exterior.coords]
        ints = [[P(c) for c in ring.coords] for ring in geom.interiors]
        return Polygon(ext, ints)
    elif geom.geom_type == 'MultiPolygon':
        polys = []
        for g in geom.geoms:
            ext = [P(c) for c in g.exterior.coords]
            ints = [[P(c) for c in ring.coords] for ring in g.interiors]
            polys.append(Polygon(ext, ints))
        return unary_union(polys).buffer(0)
    else:
        raise ValueError('Somente Polygon/MultiPolygon.')

def to_wgs84(geom, crs_from):
    from pyproj import Transformer
    tr = Transformer.from_crs(crs_from, 'epsg:4326', always_xy=True)
    def P(pt): return tr.transform(pt[0], pt[1])
    ext = [P(c) for c in geom.exterior.coords]
    ints = [[P(c) for c in ring.coords] for ring in geom.interiors]
    from shapely.geometry import Polygon
    return Polygon(ext, ints)

def _bearing_to_ccw_deg(bearing_deg: float) -> float:
    ang = 90.0 - bearing_deg
    while ang <= -180.0: ang += 360.0
    while ang > 180.0: ang -= 360.0
    return ang

def make_square_5km(center_m, azimuth_deg: float):
    from shapely.geometry import box
    from shapely.affinity import rotate, translate
    L = 5000.0
    sq = box(-L/2, -L/2, L/2, L/2)
    sq = rotate(sq, azimuth_deg, origin=(0,0), use_radians=False)
    sq = translate(sq, xoff=center_m[0], yoff=center_m[1])
    return sq

# ---------- TLE parsing ----------
def load_tle_blocks(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Arquivo TLE não encontrado: {path}")
    text = open(path, 'r', encoding='utf-8', errors='ignore').read().strip()
    lines = [ln for ln in text.splitlines() if ln.strip()]
    blocks = []
    i = 0
    while i < len(lines):
        if lines[i].startswith('1 ') and i+1 < len(lines) and lines[i+1].startswith('2 '):
            name = f'UNKNOWN-{len(blocks)+1}'
            blocks.append((name, lines[i], lines[i+1])); i += 2
        elif i+2 < len(lines) and lines[i+1].startswith('1 ') and lines[i+2].startswith('2 '):
            blocks.append((lines[i], lines[i+1], lines[i+2])); i += 3
        else:
            i += 1
    return blocks

def pick_c10(blocks):
    for b in blocks:
        head = b[0].strip().upper().replace('-', ' ')
        if head == 'GHGSAT C10':
            return b
    for b in blocks:
        head = b[0].strip().upper().replace('-', ' ')
        if 'C10' in head:
            return b
    if len(blocks)==1: return blocks[0]
    raise KeyError("Não achei TLE do GHGSAT-C10 em ./tle/ghgsat.tle")

# ---------- passes ----------
def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0088
    from math import radians, sin, cos, sqrt, atan2
    dlat = radians(lat2-lat1); dlon = radians(lon2-lon1)
    a = sin(dlat/2)**2 + cos(radians(lat1))*cos(radians(lat2))*sin(dlon/2)**2
    c = 2*atan2(sqrt(a), sqrt(1-a))
    return R*c

def weekly_passes(cent_lat, cent_lon, tle_block, start_utc, end_utc, radius_km=300.0, step_seconds=60):
    ts = load.timescale()
    sat = EarthSatellite(tle_block[1], tle_block[2], tle_block[0])
    total_seconds = int((end_utc - start_utc).total_seconds())
    N = max(1, total_seconds // step_seconds + 1)
    times = ts.from_datetimes([start_utc + timedelta(seconds=i*step_seconds) for i in range(N)])
    subs = wgs84.subpoint(sat.at(times))
    lats = subs.latitude.degrees
    lons = subs.longitude.degrees
    dists = np.array([haversine_km(cent_lat, cent_lon, lat, lon) for lat,lon in zip(lats, lons)])
    close = np.where(dists <= radius_km)[0]
    if close.size == 0: return []
    groups = []
    grp = [int(close[0])]
    for idx in close[1:]:
        if idx - grp[-1] <= 5:
            grp.append(int(idx))
        else:
            groups.append(grp); grp=[int(idx)]
    groups.append(grp)
    passes = []
    for g in groups:
        garr = np.array(g)
        local = dists[garr]
        k = garr[np.argmin(local)]
        t_mid = start_utc + timedelta(seconds=int(k*step_seconds))
        ts0 = ts.from_datetime(t_mid.replace(tzinfo=timezone.utc))
        ts1 = ts.from_datetime((t_mid + timedelta(seconds=1)).replace(tzinfo=timezone.utc))
        sp0 = wgs84.subpoint(sat.at(ts0))
        sp1 = wgs84.subpoint(sat.at(ts1))
        import math
        lat1 = math.radians(sp0.latitude.degrees); lon1 = math.radians(sp0.longitude.degrees)
        lat2 = math.radians(sp1.latitude.degrees); lon2 = math.radians(sp1.longitude.degrees)
        dlon = lon2 - lon1
        x = math.sin(dlon) * math.cos(lat2)
        y = math.cos(lat1)*math.sin(lat2) - math.sin(lat1)*math.cos(lat2)*math.cos(dlon)
        bearing = (math.degrees(math.atan2(x, y)) + 360.0) % 360.0
        kind = 'ASC' if (sp1.latitude.degrees - sp0.latitude.degrees) > 0 else 'DESC'
        passes.append({'time_utc': t_mid.isoformat(), 'bearing': bearing, 'kind': kind, 'dist_km': float(dists[k])})
    return passes

# ---------- UI (two columns only) ----------
st.title("🛰️ GHGSat C10 — Footprints semanais (ASC/DESC) — TLE offline")
st.caption("Desenhe a AOI, confirme, e veja os footprints 5×5 km de todas as passagens na próxima semana (UTC).")

col_left, col_right = st.columns([1.2, 1.2])

if 'aoi' not in st.session_state: st.session_state['aoi'] = None
if 'foot_fc' not in st.session_state: st.session_state['foot_fc'] = None

with col_left:
    st.subheader("1) Desenhar AOI")
    m1 = folium.Map(location=[-15.78, -47.93], zoom_start=4, tiles='OpenStreetMap')
    from folium.plugins import Draw
    Draw(export=False, filename='aoi.geojson',
         draw_options={'polygon': True, 'rectangle': True, 'circle': False, 'polyline': False, 'marker': False, 'circlemarker': False},
         edit_options={'edit': True, 'remove': True}).add_to(m1)
    left_map = st_folium(m1, height=560, returned_objects=['last_drawn_feature','all_drawings'])

    if st.button("Confirmar AOI", type="primary"):
        feat = None
        if left_map and left_map.get('last_drawn_feature'):
            feat = left_map['last_drawn_feature']
        elif left_map and left_map.get('all_drawings'):
            # une todos os polígonos
            feats = left_map['all_drawings'].get('features', [])
            polys = []
            for f in feats:
                g = f.get('geometry')
                if not g: continue
                t = g.get('type')
                if t in ('Polygon','MultiPolygon'):
                    polys.append(shape(g))
                elif t == 'GeometryCollection':
                    for gg in g.get('geometries', []):
                        if gg.get('type') in ('Polygon','MultiPolygon'):
                            polys.append(shape(gg))
            if polys:
                merged = unary_union(polys).buffer(0)
                feat = {'type':'Feature','geometry':mapping(merged),'properties':{}}
        if not feat:
            st.warning("Desenhe um polígono antes de confirmar.")
        else:
            st.session_state['aoi'] = feat['geometry']
            # --- compute passes & footprints ---
            try:
                blocks = load_tle_blocks(TLE_PATH)
                tle_block = pick_c10(blocks)
            except Exception as e:
                st.error(f"Erro TLE: {e}")
                st.stop()
            now = datetime.utcnow().replace(tzinfo=timezone.utc)
            end = now + timedelta(days=7)
            cent = shape(st.session_state['aoi']).centroid
            passes = weekly_passes(cent.y, cent.x, tle_block, now, end, radius_km=300.0, step_seconds=60)

            # footprints
            crs_utm = utm_crs_from_lonlat(cent.x, cent.y)
            aoi_m = to_proj(st.session_state['aoi'], crs_utm)
            c = aoi_m.centroid
            features = []
            for p in passes:
                rot = _bearing_to_ccw_deg(p['bearing'])
                foot_m = make_square_5km((c.x, c.y), rot)
                foot = to_wgs84(foot_m, crs_utm)
                features.append({'type':'Feature','geometry':mapping(foot),
                                 'properties':{'time_utc':p['time_utc'], 'kind':p['kind'], 'bearing':round(p['bearing'],2), 'dist_km':round(p['dist_km'],1)}})
            st.session_state['foot_fc'] = {'type':'FeatureCollection','features':features}
            st.success(f"Passagens encontradas: {len(features)}")

with col_right:
    st.subheader("2) Footprints semanais")
    m2 = folium.Map(location=[-15.78, -47.93], zoom_start=4, tiles='OpenStreetMap')

    if st.session_state['aoi']:
        folium.GeoJson({'type':'FeatureCollection','features':[{'type':'Feature','geometry':st.session_state['aoi']}]},
                       style_function=lambda x: {'fillColor':'#ffffff','color':'#000','weight':2,'fillOpacity':0.25}).add_to(m2)

    if st.session_state['foot_fc']:
        show_asc = st.checkbox("Mostrar ASC", value=True)
        show_desc = st.checkbox("Mostrar DESC", value=True)
        feats = [f for f in st.session_state['foot_fc']['features']
                 if (f['properties']['kind']=='ASC' and show_asc) or (f['properties']['kind']=='DESC' and show_desc)]
        fc = {'type':'FeatureCollection','features':feats}
        folium.GeoJson(fc,
                       style_function=lambda x: {'fillColor':'#88c','color':'#224','weight':2,'opacity':0.9,'fillOpacity':0.35},
                       tooltip=folium.GeoJsonTooltip(fields=['time_utc','kind','bearing','dist_km'],
                                                     aliases=['UTC','Tipo','Azimute (°)','Dist. (km)'])).add_to(m2)
        st.download_button("Baixar GeoJSON (semana)", data=json.dumps(fc), file_name="ghgsat_c10_week_footprints.geojson",
                           mime="application/geo+json")
        rows = sorted([f['properties'] for f in feats], key=lambda d: d['time_utc'])
        st.dataframe(rows, use_container_width=True)

    st_folium(m2, height=560)
