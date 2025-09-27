# -*- coding: utf-8 -*-
# app.py — GHGSat 5×5 km Planner (v8)
# - Apenas GHGSat (footprint 5 km × 5 km, 1 tile)
# - Azimute pode ser calculado automaticamente a partir do TLE mais recente (Celestrak)
# - Usuário pode escolher um satélite GHGSat na lista, ou colar TLE manualmente
# - Fallback automático para azimute manual

import json
from datetime import datetime, timezone, timedelta
from urllib.parse import urlencode

import streamlit as st
import folium
from streamlit_folium import st_folium
from shapely.geometry import Polygon, shape, mapping, box
from shapely.ops import unary_union
from shapely.affinity import rotate, translate
from pyproj import CRS, Transformer

st.set_page_config(page_title='GHGSat 5×5 Planner (TLE automático)', layout='wide')

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

def compute_track_azimuth_from_tle(tle_text: str, when_utc: datetime) -> float:
    from skyfield.api import EarthSatellite, load, wgs84
    # Parse TLE: 2 ou 3 linhas
    lines = [ln.strip() for ln in tle_text.strip().splitlines() if ln.strip()]
    if len(lines) == 2:
        name = "GHGSat"
        l1, l2 = lines
    elif len(lines) >= 3:
        name, l1, l2 = lines[0], lines[1], lines[2]
    else:
        raise ValueError("TLE inválido (forneça 2 ou 3 linhas).")
    sat = EarthSatellite(l1, l2, name)
    ts = load.timescale()
    t0 = ts.from_datetime(when_utc.replace(tzinfo=timezone.utc))

    geocentric0 = sat.at(t0)
    sp0 = wgs84.subpoint(geocentric0)
    geocentric1 = sat.at(ts.from_datetime((when_utc + timedelta(seconds=1)).replace(tzinfo=timezone.utc)))
    sp1 = wgs84.subpoint(geocentric1)

    import math
    lat1 = math.radians(sp0.latitude.degrees)
    lon1 = math.radians(sp0.longitude.degrees)
    lat2 = math.radians(sp1.latitude.degrees)
    lon2 = math.radians(sp1.longitude.degrees)
    dlon = lon2 - lon1
    x = math.sin(dlon) * math.cos(lat2)
    y = math.cos(lat1)*math.sin(lat2) - math.sin(lat1)*math.cos(lat2)*math.cos(dlon)
    bearing = (math.degrees(math.atan2(x, y)) + 360.0) % 360.0
    return bearing

def make_square_5km(center_m, azimuth_deg: float):
    L = 5000.0
    sq = box(-L/2, -L/2, L/2, L/2)
    sq = rotate(sq, azimuth_deg, origin=(0,0), use_radians=False)
    sq = translate(sq, xoff=center_m[0], yoff=center_m[1])
    return sq

# ---------- TLE fetch (Celestrak) ----------
GHGSAT_NAMES = [
    "GHGSAT-C1","GHGSAT-C2","GHGSAT-C3","GHGSAT-C4","GHGSAT-C5","GHGSAT-C6",
    "GHGSAT-C7","GHGSAT-C8","GHGSAT-C9","GHGSAT-C10","GHGSAT-C11","GHGSAT-C12",
]

def fetch_tle_from_celestrak(name_or_id: str) -> str:
    """
    Tenta baixar TLE do Celestrak usando NAME (mais comum) ou NORAD id numérico.
    URLs usadas (retornam texto TLE ou vazio):
      - https://celestrak.org/NORAD/elements/gp.php?{NAME=<>,FORMAT=TLE}
      - https://celestrak.org/NORAD/elements/gp.php?{CATNR=<num>,FORMAT=TLE}
    """
    import requests
    base = "https://celestrak.org/NORAD/elements/gp.php"
    params = {"FORMAT":"TLE"}
    # Detecta se é número (NORAD)
    is_num = name_or_id.strip().isdigit()
    if is_num:
        params["CATNR"] = name_or_id.strip()
    else:
        params["NAME"] = name_or_id.strip()
    resp = requests.get(base, params=params, timeout=10)
    if resp.status_code == 200 and resp.text and "1 " in resp.text and "2 " in resp.text:
        return resp.text.strip()
    # fallback: tenta NAME em upper
    if not is_num:
        params["NAME"] = name_or_id.strip().upper()
        resp = requests.get(base, params=params, timeout=10)
        if resp.status_code == 200 and resp.text and "1 " in resp.text and "2 " in resp.text:
            return resp.text.strip()
    raise RuntimeError(f"Não consegui TLE para '{name_or_id}'.")

# ------------- UI -------------
st.title("🛰️ GHGSat — Footprint 5 × 5 km (TLE automático)")
st.caption("Escolha um satélite GHGSat para buscar TLE no Celestrak, ou cole um TLE manualmente. O quadrado 5×5 km é rotacionado pelo azimute obtido.")

col_left, col_mid, col_right = st.columns([1.1, 0.6, 1.1])

if 'aoi' not in st.session_state: st.session_state['aoi'] = None
if 'foot' not in st.session_state: st.session_state['foot'] = None
if 'az' not in st.session_state: st.session_state['az'] = 0.0

# ---------------- AOI ----------------
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

# ---------------- Orientation ----------------
with col_mid:
    st.subheader("2) Orientação pela órbita")
    when = st.text_input("Data e hora UTC (YYYY-MM-DD HH:MM:SS)", value=datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    mode = st.radio("Como obter o azimute?", ["TLE automático (Celestrak)", "TLE manual (colar abaixo)", "Azimute manual (graus)"], index=0)

    ghg_name = st.selectbox("Satélite GHGSat", GHGSAT_NAMES, index=9)  # C10 padrão
    tle_manual = st.text_area("TLE (somente se você escolher 'TLE manual')", height=120, value="")
    az_manual = st.number_input("Azimute manual (°)", min_value=-180.0, max_value=360.0, value=0.0, step=5.0)

    plan = st.button("Planejar ➜", use_container_width=True)

# ---------------- Map out ----------------
with col_right:
    st.subheader("3) Footprint 5×5 km")
    m2 = folium.Map(location=[-15.78, -47.93], zoom_start=6, tiles='OpenStreetMap')

    if plan:
        if not st.session_state['aoi']:
            st.warning("Confirme uma AOI primeiro.")
        else:
            # Projeta AOI em metros e calcula centróide
            aoi_geo = st.session_state['aoi']
            cent_ll = shape(aoi_geo).centroid  # WGS84
            crs_utm = utm_crs_from_lonlat(cent_ll.x, cent_ll.y)
            aoi_m = to_proj(aoi_geo, crs_utm)
            c = aoi_m.centroid

            # Resolve azimute
            az = None
            try:
                when_dt = datetime.strptime(when, "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc)
            except Exception:
                when_dt = datetime.utcnow().replace(tzinfo=timezone.utc)
                st.warning("Formato de data/hora inválido. Usando agora (UTC).")

            if mode == "TLE automático (Celestrak)":
                try:
                    tle_text = fetch_tle_from_celestrak(ghg_name)
                    az = compute_track_azimuth_from_tle(tle_text, when_dt)
                except Exception as e:
                    st.error(f"Falha TLE automático: {e}")
            elif mode == "TLE manual (colar abaixo)":
                if tle_manual.strip():
                    try:
                        az = compute_track_azimuth_from_tle(tle_manual, when_dt)
                    except Exception as e:
                        st.error(f"TLE manual inválido: {e}")
                else:
                    st.warning("Cole um TLE válido.")

            if az is None and mode == "Azimute manual (graus)":
                az = float(az_manual)
            if az is None:
                st.warning("Não consegui obter azimute; usando 0°.")
                az = 0.0
            st.session_state['az'] = az

            # Footprint 5x5
            foot_m = make_square_5km((c.x, c.y), azimuth_deg=az)
            foot_wgs = to_wgs84(foot_m, crs_utm)
            st.session_state['foot'] = mapping(foot_wgs)

    # Desenha AOI
    if st.session_state['aoi']:
        folium.GeoJson({'type':'FeatureCollection','features':[{'type':'Feature','geometry':st.session_state['aoi']}]},
                       style_function=lambda x: {'fillColor':'#ffffff','color':'#000','weight':2,'fillOpacity':0.25}).add_to(m2)

    # Desenha footprint
    if st.session_state['foot']:
        folium.GeoJson({'type':'FeatureCollection','features':[{'type':'Feature','geometry':st.session_state['foot']}]} ,
                       style_function=lambda x: {'fillColor':'#6d8b9b','color':'#2b3a42','weight':2,'opacity':0.95,'fillOpacity':0.55}).add_to(m2)
        st.write(f"Azimute usado: {st.session_state['az']:.2f}°")
        st.download_button("Baixar GeoJSON (5x5 km)", data=json.dumps({'type':'FeatureCollection','features':[{'type':'Feature','geometry':st.session_state['foot']}] }),
                           file_name="ghgsat_footprint_5x5.geojson", mime="application/geo+json")

    st_folium(m2, height=560)
