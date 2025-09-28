# -*- coding: utf-8 -*-
# app.py — GHGSat 5×5 km Planner (v9c — TLE offline no repositório, sem painel de cache)
# - Lê TLE exclusivamente de ./tle/ghgsat.tle (você mantém via commit/push)
# - Remove qualquer exibição de "cache local" na UI

import os, io, re, json
from datetime import datetime, timezone, timedelta

import streamlit as st
import folium
from streamlit_folium import st_folium
from shapely.geometry import Polygon, shape, mapping, box
from shapely.ops import unary_union
from shapely.affinity import rotate, translate
from pyproj import CRS, Transformer

st.set_page_config(page_title='GHGSat 5×5 Planner (TLE offline)', layout='wide')

TLE_PATH = os.path.join(os.path.dirname(__file__), 'tle', 'ghgsat.tle')

GHGSAT_NAMES = [
    "GHGSAT-C1","GHGSAT-C2","GHGSAT-C3","GHGSAT-C4","GHGSAT-C5","GHGSAT-C6",
    "GHGSAT-C7","GHGSAT-C8","GHGSAT-C9","GHGSAT-C10","GHGSAT-C11","GHGSAT-C12",
]

# ------------- Helpers -------------
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

# ------------- TLE parsing (offline) -------------
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

def index_tles(blocks):
    by_name, by_norad = {}, {}
    for name, l1, l2 in blocks:
        name_key = name.strip().upper().replace('-', ' ')
        by_name[name_key] = (name, l1, l2)
        m = re.match(r'^1\s+(\d+)\s', l1.strip())
        if m:
            by_norad[m.group(1)] = (name, l1, l2)
    return by_name, by_norad

def get_tle_from_cache(name_or_catnr: str):
    blocks = load_tle_file(TLE_PATH)
    by_name, by_norad = index_tles(blocks)
    q_raw = name_or_catnr.strip()
    q = q_raw.upper()
    q_norm = q.replace('-', ' ').replace('  ', ' ')
    if q in by_name: return by_name[q]
    if q_norm in by_name: return by_name[q_norm]
    import re as _re
    m = _re.search(r'(C\s*\d{1,2})', q_norm)
    if m:
        ctag = m.group(1).replace(' ', '')
        for key in list(by_name.keys()):
            key_norm = key.replace('-', ' ').replace('  ', ' ')
            if ctag in key_norm.replace(' ', ''):
                return by_name[key]
    for key in list(by_name.keys()):
        key_norm = key.replace('-', ' ').replace('  ', ' ')
        if q_norm in key_norm or q in key:
            return by_name[key]
    if q.isdigit() and q in by_norad:
        return by_norad[q]
    unnamed = [b for b in blocks if b[0].upper().startswith('UNKNOWN-')]
    if len(blocks) == 1:
        return blocks[0]
    if len(unnamed) == 1 and len(blocks) == 1:
        return unnamed[0]
    raise KeyError(f"Não achei TLE no arquivo para: {name_or_catnr}. "
                   f"Edite ./tle/ghgsat.tle e garanta um bloco com o nome 'GHGSAT-C10' ou informe o NORAD CATNR.")

def compute_track_azimuth_from_tle_block(tle_block, when_utc: datetime) -> float:
    from skyfield.api import EarthSatellite, load, wgs84
    name, l1, l2 = tle_block
    sat = EarthSatellite(l1, l2, name)
    ts = load.timescale()
    t0 = ts.from_datetime(when_utc.replace(tzinfo=timezone.utc))
    g0 = sat.at(t0); sp0 = wgs84.subpoint(g0)
    g1 = sat.at(ts.from_datetime((when_utc + timedelta(seconds=1)).replace(tzinfo=timezone.utc))); sp1 = wgs84.subpoint(g1)
    import math
    lat1 = math.radians(sp0.latitude.degrees); lon1 = math.radians(sp0.longitude.degrees)
    lat2 = math.radians(sp1.latitude.degrees); lon2 = math.radians(sp1.longitude.degrees)
    dlon = lon2 - lon1
    x = math.sin(dlon) * math.cos(lat2)
    y = math.cos(lat1)*math.sin(lat2) - math.sin(lat1)*math.cos(lat2)*math.cos(dlon)
    bearing = (math.degrees(math.atan2(x, y)) + 360.0) % 360.0
    return bearing

# ------------- UI -------------
st.title("🛰️ GHGSat — Footprint 5 × 5 km (TLE offline)")
st.caption("O app lê o TLE de ./tle/ghgsat.tle no repositório. Atualize esse arquivo quando quiser (commit/push).")

col_left, col_mid, col_right = st.columns([1.1, 0.7, 1.1])

if 'aoi' not in st.session_state: st.session_state['aoi'] = None
if 'foot' not in st.session_state: st.session_state['foot'] = None
if 'az' not in st.session_state: st.session_state['az'] = 0.0
if 'az_is_compass' not in st.session_state: st.session_state['az_is_compass'] = True

# Esquerda: AOI
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
            st.session_state['foot'] = None
            st.success("AOI confirmada!")

# Meio: orientação
with col_mid:
    st.subheader("2) Orientação pela órbita (via TLE local)")
    when = st.text_input("Data e hora UTC (YYYY-MM-DD HH:MM:SS)", value=datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    mode = st.radio("Como obter o azimute?", ["Do arquivo TLE local", "TLE colado (manual)", "Azimute manual (geométrico)"], index=0)

    ghg_name = st.selectbox("Satélite GHGSat", GHGSAT_NAMES, index=9)  # C10 padrão
    norad_id_opt = st.text_input("NORAD CATNR (opcional)", value="")
    tle_manual = st.text_area("TLE (somente se escolher 'TLE colado')", height=120, value="")
    az_manual = st.number_input("Azimute manual (°) — 0=E, CCW", min_value=-180.0, max_value=360.0, value=0.0, step=5.0)

    plan = st.button("Planejar ➜", use_container_width=True)

# Direita: saída
with col_right:
    st.subheader("3) Footprint 5×5 km")
    m2 = folium.Map(location=[-15.78, -47.93], zoom_start=6, tiles='OpenStreetMap')

    if plan:
        if not st.session_state['aoi']:
            st.warning("Confirme uma AOI primeiro.")
        else:
            aoi_geo = st.session_state['aoi']
            cent_ll = shape(aoi_geo).centroid
            crs_utm = utm_crs_from_lonlat(cent_ll.x, cent_ll.y)
            aoi_m = to_proj(aoi_geo, crs_utm)
            c = aoi_m.centroid

            try:
                when_dt = datetime.strptime(when, "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc)
            except Exception:
                st.warning("Data/hora inválida; usando agora (UTC).")
                when_dt = datetime.utcnow().replace(tzinfo=timezone.utc)

            az = None; az_is_compass = False

            if mode == "Do arquivo TLE local":
                try:
                    target = (norad_id_opt.strip() or ghg_name).strip()
                    tle_block = get_tle_from_cache(target)
                    az = compute_track_azimuth_from_tle_block(tle_block, when_dt)  # bearing bússola
                    az_is_compass = True
                except Exception as e:
                    st.error(f"Falha ao ler/calcular com TLE local: {e}")
                    st.info("Dica: edite ./tle/ghgsat.tle e garanta um bloco com o nome 'GHGSAT-C10' na linha 1, seguido das linhas '1 ...' e '2 ...'. Ou preencha o NORAD CATNR.")

            elif mode == "TLE colado (manual)":
                if tle_manual.strip():
                    try:
                        lines = [ln for ln in tle_manual.strip().splitlines() if ln.strip()]
                        if len(lines) >= 3:
                            tle_block = (lines[0], lines[1], lines[2])
                        elif len(lines) == 2:
                            tle_block = ("GHGSAT", lines[0], lines[1])
                        else:
                            raise ValueError("Forneça 2 ou 3 linhas de TLE.")
                        az = compute_track_azimuth_from_tle_block(tle_block, when_dt)
                        az_is_compass = True
                    except Exception as e:
                        st.error(f"TLE manual inválido: {e}")
                else:
                    st.warning("Cole um TLE válido.")

            elif mode == "Azimute manual (geométrico)":
                az = float(az_manual); az_is_compass = False

            if az is None:
                st.warning("Não consegui obter azimute; usando 0° (Leste).")
                az = 0.0

            rot_deg = _bearing_to_ccw_deg(az) if az_is_compass else az
            foot_m = make_square_5km((c.x, c.y), azimuth_deg=rot_deg)
            foot_wgs = to_wgs84(foot_m, crs_utm)
            st.session_state['foot'] = mapping(foot_wgs)
            st.session_state['az'] = az
            st.session_state['az_is_compass'] = az_is_compass

    if st.session_state['aoi']:
        folium.GeoJson({'type':'FeatureCollection','features':[{'type':'Feature','geometry':st.session_state['aoi']}]},
                       style_function=lambda x: {'fillColor':'#ffffff','color':'#000','weight':2,'fillOpacity':0.25}).add_to(m2)

    if st.session_state['foot']:
        folium.GeoJson({'type':'FeatureCollection','features':[{'type':'Feature','geometry':st.session_state['foot']}]} ,
                       style_function=lambda x: {'fillColor':'#6d8b9b','color':'#2b3a42','weight':2,'opacity':0.95,'fillOpacity':0.55}).add_to(m2)
        st.write(f"Azimute usado: {st.session_state['az']:.2f}° " + ("(bússola/TLE)" if st.session_state['az_is_compass'] else "(geométrico CCW)"))
        st.download_button("Baixar GeoJSON (5x5 km)", data=json.dumps({'type':'FeatureCollection','features':[{'type':'Feature','geometry':st.session_state['foot']}] }),
                           file_name="ghgsat_footprint_5x5.geojson", mime="application/geo+json")

    st_folium(m2, height=560)
