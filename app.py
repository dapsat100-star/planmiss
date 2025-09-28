# -*- coding: utf-8 -*-
# GHGSat C10 — Footprints semanais (v9e-fix-ui+colors, 2 colunas, TLE offline, UI grande)

import os, json, math
from datetime import datetime, timezone, timedelta

import streamlit as st
import folium
from streamlit_folium import st_folium
from shapely.geometry import shape, mapping, box, Polygon
from shapely.ops import unary_union
from shapely.affinity import rotate, translate
from pyproj import CRS, Transformer
from skyfield.api import EarthSatellite, load, wgs84
import numpy as np

# ---------- LAYOUT ----------
st.set_page_config(page_title='GHGSat C10 — Semana (TLE offline)', layout='wide')
st.markdown("""
<style>
.block-container {padding-top: 0.6rem; padding-bottom: 0.4rem; max-width: 95vw;}
.mapbox { position: relative; }
.mapbadge {
  position: absolute; z-index: 20; top: 12px; left: 14px;
  padding: .35rem .7rem; border-radius: .6rem;
  background: rgba(15,15,20,.65); color: #fff; font-weight: 700; font-size: 0.95rem;
  box-shadow: 0 6px 18px rgba(0,0,0,.18); backdrop-filter: blur(4px);
}
.legendbadge {
  position: absolute; z-index: 20; top: 12px; right: 14px;
  padding: .35rem .6rem; border-radius: .6rem;
  background: rgba(15,15,20,.55); color: #fff; font-size: 0.85rem;
  box-shadow: 0 6px 18px rgba(0,0,0,.18); backdrop-filter: blur(4px);
}
.legendbadge .dot {display:inline-block; width:10px; height:10px; border-radius:50%; margin:0 6px 0 10px;}
.dot.asc {background:#2ecc71;}
.dot.desc {background:#e74c3c;}
.leaflet-control { filter: drop-shadow(0 2px 6px rgba(0,0,0,.25)); }
</style>
""", unsafe_allow_html=True)

TLE_PATH = os.path.join(os.path.dirname(__file__), 'tle', 'ghgsat.tle')

# --------- helpers geo ---------
def utm_crs_from_lonlat(lon: float, lat: float):
    zone = int((lon + 180) // 6) + 1
    south = lat < 0
    return CRS.from_dict({'proj': 'utm', 'zone': zone, 'south': south})

def to_proj(geom_geojson, crs_to: CRS):
    g = shape(geom_geojson); tr = Transformer.from_crs('epsg:4326', crs_to, always_xy=True)
    def P(pt): return tr.transform(pt[0], pt[1])
    if g.geom_type == 'Polygon':
        ext = [P(c) for c in g.exterior.coords]
        holes = [[P(c) for c in r.coords] for r in g.interiors]
        return Polygon(ext, holes)
    elif g.geom_type == 'MultiPolygon':
        parts = []
        for gg in g.geoms:
            ext = [P(c) for c in gg.exterior.coords]
            holes = [[P(c) for c in r.coords] for r in gg.interiors]
            parts.append(Polygon(ext, holes))
        return unary_union(parts).buffer(0)
    else:
        raise ValueError('Somente Polygon/MultiPolygon.')

def to_wgs84(geom: Polygon, crs_from: CRS):
    tr = Transformer.from_crs(crs_from, 'epsg:4326', always_xy=True)
    def P(pt): return tr.transform(pt[0], pt[1])
    ext = [P(c) for c in geom.exterior.coords]
    holes = [[P(c) for c in r.coords] for r in geom.interiors]
    return Polygon(ext, holes)

def _bearing_to_ccw_deg(brg: float) -> float:
    ang = 90.0 - brg
    while ang <= -180.0: ang += 360.0
    while ang > 180.0: ang -= 360.0
    return ang

def make_square_5km(center_m, azimuth_deg: float):
    L = 5000.0
    sq = box(-L/2, -L/2, L/2, L/2)
    sq = rotate(sq, azimuth_deg, origin=(0,0), use_radians=False)
    return translate(sq, xoff=center_m[0], yoff=center_m[1])

# --------- TLE (offline) ----------
def load_tle_blocks(path: str):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Arquivo TLE não encontrado: {path}")
    lines = [ln for ln in open(path, 'r', encoding='utf-8', errors='ignore').read().splitlines() if ln.strip()]
    blocks = []; i = 0
    while i < len(lines):
        if lines[i].startswith('1 ') and i+1 < len(lines) and lines[i+1].startswith('2 '):
            blocks.append((f'UNKNOWN-{len(blocks)+1}', lines[i], lines[i+1])); i += 2
        elif i+2 < len(lines) and lines[i+1].startswith('1 ') and lines[i+2].startswith('2 '):
            blocks.append((lines[i], lines[i+1], lines[i+2])); i += 3
        else:
            i += 1
    return blocks

def pick_c10(blocks):
    for b in blocks:
        head = b[0].strip().upper().replace('-', ' ')
        if head == 'GHGSAT C10': return b
    for b in blocks:
        if 'C10' in b[0].upper(): return b
    if len(blocks) == 1: return blocks[0]
    raise KeyError("Não achei TLE do GHGSAT-C10 em ./tle/ghgsat.tle")

# --------- passes (1 semana) ----------
def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0088
    from math import radians, sin, cos, sqrt, atan2
    dlat = radians(lat2-lat1); dlon = radians(lon2-lon1)
    a = sin(dlat/2)**2 + cos(radians(lat1))*cos(radians(lat2))*sin(dlon/2)**2
    return 2*R*atan2(sqrt(a), sqrt(1-a))

def weekly_passes(cent_lat, cent_lon, tle_block, start_utc, end_utc, radius_km=300.0, step_seconds=60):
    ts = load.timescale()
    sat = EarthSatellite(tle_block[1], tle_block[2], tle_block[0])
    total = int((end_utc - start_utc).total_seconds()); N = max(1, total // step_seconds + 1)
    times = ts.from_datetimes([start_utc + timedelta(seconds=i*step_seconds) for i in range(N)])
    subs = wgs84.subpoint(sat.at(times))
    lats = subs.latitude.degrees; lons = subs.longitude.degrees
    dists = np.array([haversine_km(cent_lat, cent_lon, la, lo) for la,lo in zip(lats, lons)])
    close = np.where(dists <= radius_km)[0]
    if close.size == 0: return []
    groups = []; grp = [int(close[0])]
    for idx in close[1:]:
        if idx - grp[-1] <= 5: grp.append(int(idx))
        else: groups.append(grp); grp=[int(idx)]
    groups.append(grp)

    passes = []
    for g in groups:
        garr = np.array(g); local = dists[garr]; k = garr[np.argmin(local)]
        t_mid = start_utc + timedelta(seconds=int(k*step_seconds))
        ts0 = ts.from_datetime(t_mid.replace(tzinfo=timezone.utc))
        ts1 = ts.from_datetime((t_mid + timedelta(seconds=1)).replace(tzinfo=timezone.utc))
        sp0 = wgs84.subpoint(sat.at(ts0)); sp1 = wgs84.subpoint(sat.at(ts1))
        lat1 = math.radians(sp0.latitude.degrees); lon1 = math.radians(sp0.longitude.degrees)
        lat2 = math.radians(sp1.latitude.degrees); lon2 = math.radians(sp1.longitude.degrees)
        dlon = lon2 - lon1
        x = math.sin(dlon) * math.cos(lat2)
        y = math.cos(lat1)*math.sin(lat2) - math.sin(lat1)*math.cos(lat2)*math.cos(dlon)
        bearing = (math.degrees(math.atan2(x, y)) + 360.0) % 360.0
        kind = 'ASC' if (sp1.latitude.degrees - sp0.latitude.degrees) > 0 else 'DESC'
        passes.append({'time_utc': t_mid.isoformat(), 'bearing': bearing, 'kind': kind, 'dist_km': float(dists[k])})
    return passes

# --------- UI (2 colunas, mapas grandes com badge/legenda) ----------
st.title("🛰️ GHGSat C10 — Footprints semanais (ASC/DESC) — TLE offline")

col_left, col_right = st.columns([1, 1], gap="medium")
MAP_HEIGHT = 820  # px

if 'aoi' not in st.session_state: st.session_state['aoi'] = None
if 'foot_fc' not in st.session_state: st.session_state['foot_fc'] = None

with col_left:
    st.markdown('<div class="mapbox"><div class="mapbadge">1) Desenhar AOI</div>', unsafe_allow_html=True)
    m1 = folium.Map(location=[-15.78, -47.93], zoom_start=4, tiles='CartoDB positron')
    from folium.plugins import Draw
    Draw(export=False, filename='aoi.geojson',
         draw_options={'polygon': True, 'rectangle': True, 'circle': False, 'polyline': False, 'marker': False, 'circlemarker': False},
         edit_options={'edit': True, 'remove': True}).add_to(m1)
    left_map = st_folium(m1, height=MAP_HEIGHT, returned_objects=['last_drawn_feature','all_drawings'])
    st.markdown('</div>', unsafe_allow_html=True)

    if st.button("✅ Confirmar AOI", use_container_width=True):
        feat = None
        if left_map and left_map.get('last_drawn_feature'):
            feat = left_map['last_drawn_feature']
        if feat is None:
            drawings = (left_map or {}).get('all_drawings')
            polys = []
            if drawings:
                if isinstance(drawings, dict) and drawings.get('type') == 'FeatureCollection':
                    features = drawings.get('features', [])
                elif isinstance(drawings, dict) and drawings.get('features'):
                    features = drawings.get('features', [])
                elif isinstance(drawings, list):
                    features = drawings
                else:
                    features = []
                for f in features:
                    g = f.get('geometry') if isinstance(f, dict) else None
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
            st.warning("Desenhe um polígono e confirme.")
        else:
            st.session_state['aoi'] = feat['geometry']
            # calcular passagens + footprints
            blocks = load_tle_blocks(TLE_PATH)
            tle_block = pick_c10(blocks)
            now = datetime.utcnow().replace(tzinfo=timezone.utc)
            end = now + timedelta(days=7)
            cent = shape(st.session_state['aoi']).centroid
            passes = weekly_passes(cent.y, cent.x, tle_block, now, end, radius_km=300.0, step_seconds=60)

            crs_utm = utm_crs_from_lonlat(cent.x, cent.y)
            aoi_m = to_proj(st.session_state['aoi'], crs_utm); c = aoi_m.centroid
            feats = []
            for p in passes:
                rot = _bearing_to_ccw_deg(p['bearing'])
                foot_m = make_square_5km((c.x, c.y), rot)
                foot = to_wgs84(foot_m, crs_utm)
                feats.append({'type':'Feature','geometry':mapping(foot),
                              'properties':{'time_utc':p['time_utc'],'kind':p['kind'],
                                            'bearing':round(p['bearing'],2),'dist_km':round(p['dist_km'],1)}})
            st.session_state['foot_fc'] = {'type':'FeatureCollection','features':feats}
            st.success(f"Passagens encontradas: {len(feats)}")

with col_right:
    st.markdown('<div class="mapbox"><div class="mapbadge">2) Footprints semanais</div>'
                '<div class="legendbadge"><span class="dot asc"></span>ASC  '
                '<span class="dot desc"></span>DESC</div>', unsafe_allow_html=True)

    m2 = folium.Map(location=[-15.78, -47.93], zoom_start=4, tiles='CartoDB positron')

    # AOI
    if st.session_state['aoi']:
        folium.GeoJson({'type':'FeatureCollection','features':[{'type':'Feature','geometry':st.session_state['aoi']}]},
                       style_function=lambda x: {'fillColor':'#ffffff','color':'#000','weight':2,'fillOpacity':0.25}).add_to(m2)

    # Footprints com cores por tipo
    if st.session_state['foot_fc']:
        def style_fn(feat):
            kind = (feat['properties'].get('kind') or '').upper()
            if kind == 'ASC':
                fill, stroke = '#2ecc71', '#1e874b'   # verde
            else:
                fill, stroke = '#e74c3c', '#b03a2e'   # vermelho
            return {'fillColor': fill, 'color': stroke, 'weight': 2, 'opacity': 0.9, 'fillOpacity': 0.35}

        folium.GeoJson(st.session_state['foot_fc'],
                       style_function=style_fn,
                       tooltip=folium.GeoJsonTooltip(fields=['time_utc','kind','bearing','dist_km'],
                                                     aliases=['UTC','Tipo','Azimute (°)','Dist. (km)'])).add_to(m2)

    st_folium(m2, height=MAP_HEIGHT)
    st.markdown('</div>', unsafe_allow_html=True)

    if st.session_state.get('foot_fc'):
        st.download_button("⬇️ Baixar GeoJSON (semana)", data=json.dumps(st.session_state['foot_fc']),
                           file_name="ghgsat_c10_week_footprints.geojson", mime="application/geo+json")
        with st.expander("Ver lista de passagens (UTC)"):
            rows = sorted([f['properties'] for f in st.session_state['foot_fc']['features']], key=lambda d: d['time_utc'])
            st.dataframe(rows, use_container_width=True)
