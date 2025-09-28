# -*- coding: utf-8 -*-
# app.py — GHGSat 5×5 km Planner (v9d)
# - Somente GHGSAT-C10
# - Sem coluna de parâmetros; após confirmar a AOI, calcula TODAS as passagens em 1 semana (UTC)
# - Gera footprints 5×5 km para passagens ascendentes e descendentes a partir de TLE local (./tle/ghgsat.tle)

import os, re, json, math
from datetime import datetime, timezone, timedelta

import streamlit as st
import folium
from streamlit_folium import st_folium
from shapely.geometry import Polygon, shape, mapping, box
from shapely.ops import unary_union
from shapely.affinity import rotate, translate
from pyproj import CRS, Transformer

st.set_page_config(page_title='GHGSat C10 — Passes 1 semana (TLE offline)', layout='wide')

TLE_PATH = os.path.join(os.path.dirname(__file__), 'tle', 'ghgsat.tle')
SAT_NAME = "GHGSAT-C10"  # fixo

# ---------------- Helpers ----------------
def _extract_features(all_drawings_obj):
    feats = []
    if not all_drawings_obj:
        return feats
    if isinstance(all_drawings_obj, dict):
        if all_drawings_obj.get('type') == 'FeatureCollection' and all_drawings_obj.get('features'):
            return list(all_drawings_obj['features'])
        if all_drawings_obj.get('features'):
            return list(all_drawings_obj['features'])
    if isinstance(all_drawings_obj, list):
        for item in all_drawings_obj:
            if isinstance(item, dict) and item.get('geometry'):
                feats.append(item)
    return feats

def utm_crs_from_lonlat(lon: float, lat: float) -> CRS:
    zone = int((lon + 180) // 6) + 1
    south = lat < 0
    return CRS.from_dict({'proj': 'utm', 'zone': zone, 'south': south})

def to_proj(geom_geojson, crs_to: CRS):
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
        raise ValueError('Somente Polygon/MultiPolygon são suportados.')

def to_wgs84(geom: Polygon, crs_from: CRS):
    tr = Transformer.from_crs(crs_from, 'epsg:4326', always_xy=True)
    def P(pt): return tr.transform(pt[0], pt[1])
    ext = [P(c) for c in geom.exterior.coords]
    ints = [[P(c) for c in ring.coords] for ring in geom.interiors]
    return Polygon(ext, ints)

def _bearing_to_ccw_deg(bearing_deg: float) -> float:
    ang = 90.0 - bearing_deg
    while ang <= -180.0: ang += 360.0
    while ang > 180.0: ang -= 360.0
    return ang

def make_square_5km(center_m, azimuth_deg: float):
    L = 5000.0
    sq = box(-L/2, -L/2, L/2, L/2)
    sq = rotate(sq, azimuth_deg, origin=(0,0), use_radians=False)
    sq = translate(sq, xoff=center_m[0], yoff=center_m[1])
    return sq

# ---------------- TLE utils (offline) ----------------
def load_tle_file(path: str):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Arquivo TLE não encontrado: {path}")
    text = open(path, 'r', encoding='utf-8', errors='ignore').read().strip()
    lines = [ln for ln in text.splitlines() if ln.strip()]
    blocks = []
    i = 0
    while i < len(lines):
        if lines[i].startswith('1 ') and i+1 < len(lines) and lines[i+1].startswith('2 '):
            name = f'UNKNOWN-{len(blocks)+1}'
            l1, l2 = lines[i], lines[i+1]
            blocks.append((name, l1, l2)); i += 2
        elif i+2 < len(lines) and lines[i+1].startswith('1 ') and lines[i+2].startswith('2 '):
            name, l1, l2 = lines[i], lines[i+1], lines[i+2]
            blocks.append((name, l1, l2)); i += 3
        else:
            i += 1
    return blocks

def pick_tle_for_c10(blocks):
    # tenta cabeça igual; senão aceita substring "C10"
    for name,l1,l2 in blocks:
        head = name.strip().upper().replace('-', ' ')
        if head == 'GHGSAT C10':
            return (name,l1,l2)
    for name,l1,l2 in blocks:
        head = name.strip().upper().replace('-', ' ')
        if 'C10' in head:
            return (name,l1,l2)
    # fallback: se só tiver um bloco, usa
    if len(blocks)==1:
        return blocks[0]
    raise KeyError("Não achei TLE do GHGSAT-C10 em ./tle/ghgsat.tle")

# ---------------- Pass planner ----------------
def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0088
    from math import radians, sin, cos, sqrt, atan2
    dlat = radians(lat2-lat1)
    dlon = radians(lon2-lon1)
    a = sin(dlat/2)**2 + cos(radians(lat1))*cos(radians(lat2))*sin(dlon/2)**2
    c = 2*atan2(sqrt(a), sqrt(1-a))
    return R*c

def find_passes_week(cent_lat, cent_lon, tle_block, start_utc: datetime, end_utc: datetime,
                     radius_km=250.0, step_seconds=60):
    """Varre em passos de 60s; encontra janelas em que o subponto passa perto do alvo.
       radius_km controla o 'quão perto' do centróide; 250 km é folgado pra garantir detecção.
       Retorna lista de dicts: {time, bearing, kind('ASC'/'DESC'), dist_km}
    """
    from skyfield.api import EarthSatellite, load, wgs84
    ts = load.timescale()
    name,l1,l2 = tle_block
    sat = EarthSatellite(l1,l2,name)
    # times array
    total_seconds = int((end_utc - start_utc).total_seconds())
    N = max(1, total_seconds // step_seconds + 1)
    times = ts.from_datetimes([start_utc + timedelta(seconds=i*step_seconds) for i in range(N)])
    geos = sat.at(times)
    subs = wgs84.subpoint(geos)
    lats = subs.latitude.degrees
    lons = subs.longitude.degrees
    # compute distances
    import numpy as np
    dists = np.array([haversine_km(cent_lat, cent_lon, lat, lon) for lat,lon in zip(lats, lons)])
    close = np.where(dists <= radius_km)[0]
    passes = []
    if close.size == 0:
        return passes
    # group contiguous indices into clusters
    groups = []
    grp = [int(close[0])]
    for idx in close[1:]:
        if idx - grp[-1] <= 5:  # <=5 min entre pontos pertence ao mesmo passe
            grp.append(int(idx))
        else:
            groups.append(grp); grp=[int(idx)]
    groups.append(grp)
    for g in groups:
        # pick min distance index in this group
        import numpy as np
        garr = np.array(g)
        local = dists[garr]
        k = garr[np.argmin(local)]
        t_mid = start_utc + timedelta(seconds=int(k*step_seconds))
        # bearing: use t_mid and t_mid+1s
        ts0 = ts.from_datetime(t_mid.replace(tzinfo=timezone.utc))
        ts1 = ts.from_datetime((t_mid + timedelta(seconds=1)).replace(tzinfo=timezone.utc))
        sp0 = wgs84.subpoint(sat.at(ts0))
        sp1 = wgs84.subpoint(sat.at(ts1))
        # compute bearing/bússola
        lat1 = math.radians(sp0.latitude.degrees); lon1 = math.radians(sp0.longitude.degrees)
        lat2 = math.radians(sp1.latitude.degrees); lon2 = math.radians(sp1.longitude.degrees)
        dlon = lon2 - lon1
        x = math.sin(dlon) * math.cos(lat2)
        y = math.cos(lat1)*math.sin(lat2) - math.sin(lat1)*math.cos(lat2)*math.cos(dlon)
        bearing = (math.degrees(math.atan2(x, y)) + 360.0) % 360.0
        # ascending/descending by latitude trend
        kind = 'ASC' if (sp1.latitude.degrees - sp0.latitude.degrees) > 0 else 'DESC'
        passes.append({
            'time_utc': t_mid.isoformat(),
            'bearing': bearing,
            'kind': kind,
            'dist_km': float(dists[k])
        })
    return passes

# ---------------- UI ----------------
st.title("🛰️ GHGSat C10 — Passagens (1 semana) com footprints 5×5 km — TLE offline")
st.caption("Desenhe e confirme a AOI. O app usa GHGSAT-C10 e calcula as passagens (ASC/DESC) para 7 dias a partir de agora (UTC).")

col_left, col_right = st.columns([1.2, 1.2])

if 'aoi' not in st.session_state: st.session_state['aoi'] = None
if 'passes' not in st.session_state: st.session_state['passes'] = []
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
        elif left_map and _extract_features(left_map.get('all_drawings')):
            feats = _extract_features(left_map.get('all_drawings'))
            polys = []
            for f in feats:
                g = f.get('geometry'); 
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
            # compute passes immediately
            aoi_geo = st.session_state['aoi']
            cent_ll = shape(aoi_geo).centroid
            # load TLE and pick C10
            blocks = load_tle_file(TLE_PATH)
            try:
                tle_block = pick_tle_for_c10(blocks)
            except Exception as e:
                st.error(f"Falha ao ler TLE C10 em ./tle/ghgsat.tle: {e}")
                st.stop()
            now = datetime.utcnow().replace(tzinfo=timezone.utc)
            end = now + timedelta(days=7)
            passes = find_passes_week(cent_ll.y, cent_ll.x, tle_block, now, end,
                                      radius_km=300.0, step_seconds=60)
            st.session_state['passes'] = passes
            # build footprints FeatureCollection
            crs_utm = utm_crs_from_lonlat(cent_ll.x, cent_ll.y)
            aoi_m = to_proj(aoi_geo, crs_utm); c = aoi_m.centroid
            feats = []
            for p in passes:
                rot = _bearing_to_ccw_deg(p['bearing'])
                foot_m = make_square_5km((c.x, c.y), rot)
                foot = to_wgs84(foot_m, crs_utm)
                feats.append({'type':'Feature','geometry':mapping(foot),
                              'properties':{'time_utc':p['time_utc'], 'kind':p['kind'], 'bearing':round(p['bearing'],2), 'dist_km':round(p['dist_km'],1)}})
            st.session_state['foot_fc'] = {'type':'FeatureCollection','features':feats}
            st.success(f"Encontradas {len(passes)} passagens na próxima semana.")

with col_right:
    st.subheader("2) Footprints (ASC/DESC) na semana")
    m2 = folium.Map(location=[-15.78, -47.93], zoom_start=4, tiles='OpenStreetMap')

    # Desenha AOI
    if st.session_state['aoi']:
        folium.GeoJson({'type':'FeatureCollection','features':[{'type':'Feature','geometry':st.session_state['aoi']}]},
                       style_function=lambda x: {'fillColor':'#ffffff','color':'#000','weight':2,'fillOpacity':0.25}).add_to(m2)

    if st.session_state['foot_fc']:
        show_asc = st.checkbox("Mostrar passagens ascendentes (ASC)", value=True)
        show_desc = st.checkbox("Mostrar passagens descendentes (DESC)", value=True)
        feats = []
        for f in st.session_state['foot_fc']['features']:
            if (f['properties']['kind']=='ASC' and show_asc) or (f['properties']['kind']=='DESC' and show_desc):
                feats.append(f)
        fc = {'type':'FeatureCollection','features':feats}
        folium.GeoJson(fc,
                       name='Footprints semana',
                       style_function=lambda x: {'fillColor':'#88c','color':'#224','weight':2,'opacity':0.9,'fillOpacity':0.35},
                       tooltip=folium.GeoJsonTooltip(fields=['time_utc','kind','bearing','dist_km'],
                                                     aliases=['UTC','Tipo','Azimute (°)','Dist. ao alvo (km)'])).add_to(m2)
        st.download_button("Baixar GeoJSON (semana)", data=json.dumps(fc), file_name="ghgsat_c10_week_footprints.geojson",
                           mime="application/geo+json")

        # tabela simples
        rows = sorted([f['properties'] for f in st.session_state['foot_fc']['features']], key=lambda d: d['time_utc'])
        st.write(f"Total de passagens: {len(rows)}")
        st.dataframe(rows, use_container_width=True)

    st_folium(m2, height=560)
