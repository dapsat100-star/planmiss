# -*- coding: utf-8 -*-
# app.py — GHGSat 5x5 km Footprint Planner (Streamlit + Folium)

import json
from typing import List, Tuple
from dataclasses import dataclass

import streamlit as st
import folium
from streamlit_folium import st_folium
from shapely.geometry import Polygon, shape, box, mapping
from shapely.ops import unary_union
from pyproj import CRS, Transformer

# ==========================
# Helper para extrair feições desenhadas
# ==========================
def _extract_features(all_drawings_obj):
    feats = []
    if not all_drawings_obj:
        return feats
    # Caso 1: FeatureCollection
    if isinstance(all_drawings_obj, dict):
        if all_drawings_obj.get("type") == "FeatureCollection" and all_drawings_obj.get("features"):
            return list(all_drawings_obj["features"])
        if all_drawings_obj.get("features"):
            return list(all_drawings_obj["features"])
    # Caso 2: lista de features
    if isinstance(all_drawings_obj, list):
        for item in all_drawings_obj:
            if isinstance(item, dict) and item.get("geometry"):
                feats.append(item)
    return feats

# ==========================
# Configuração inicial
# ==========================
st.set_page_config(page_title="Planejador GHGSat 5x5 km", layout="wide")

# ==========================
# Funções de projeção
# ==========================
def utm_crs_from_lonlat(lon: float, lat: float) -> CRS:
    zone = int((lon + 180) // 6) + 1
    south = lat < 0
    return CRS.from_dict({"proj": "utm", "zone": zone, "south": south})

def to_proj(geom_geojson, crs_to: CRS):
    geom = shape(geom_geojson)
    transformer = Transformer.from_crs("epsg:4326", crs_to, always_xy=True)
    def _reproj_coords(coords):
        return transformer.transform(coords[0], coords[1])
    if geom.geom_type == "Polygon":
        exterior = [_reproj_coords(c) for c in geom.exterior.coords]
        interiors = [[_reproj_coords(c) for c in ring.coords] for ring in geom.interiors]
        return Polygon(exterior, interiors)
    elif geom.geom_type == "MultiPolygon":
        polys = []
        for g in geom.geoms:
            exterior = [_reproj_coords(c) for c in g.exterior.coords]
            interiors = [[_reproj_coords(c) for c in ring.coords] for ring in g.interiors]
            polys.append(Polygon(exterior, interiors))
        return unary_union(polys).buffer(0)
    else:
        raise ValueError("Somente Polygon/MultiPolygon são suportados.")

def to_wgs84(geom: Polygon, crs_from: CRS):
    transformer = Transformer.from_crs(crs_from, "epsg:4326", always_xy=True)
    def _reproj_coords(coords):
        return transformer.transform(coords[0], coords[1])
    exterior = [_reproj_coords(c) for c in geom.exterior.coords]
    interiors = [[_reproj_coords(c) for c in ring.coords] for ring in geom.interiors]
    return Polygon(exterior, interiors)

# ==========================
# Função para gerar grid 5x5 km
# ==========================
@dataclass
class GridResult:
    tiles_geojson: dict
    n_tiles: int
    total_area_km2: float

def generate_5km_grid_covering(aoi_geojson: dict) -> GridResult:
    centroid = shape(aoi_geojson).centroid
    crs_utm = utm_crs_from_lonlat(centroid.x, centroid.y)
    aoi_m = to_proj(aoi_geojson, crs_utm)

    cell = 5000.0  # metros
    minx, miny, maxx, maxy = aoi_m.bounds
    minx -= cell; miny -= cell; maxx += cell; maxy += cell

    def snap(v, size):
        return (int(v // size)) * size

    start_x = snap(minx, cell)
    start_y = snap(miny, cell)

    tiles = []
    y = start_y
    while y < maxy:
        x = start_x
        while x < maxx:
            tile = box(x, y, x+cell, y+cell)
            if tile.intersects(aoi_m):
                tiles.append(tile.intersection(aoi_m.buffer(0)))
            x += cell
        y += cell

    feats = []
    total_area = 0.0
    for t in tiles:
        if t.is_empty:
            continue
        t_wgs = to_wgs84(t, crs_utm)
        feats.append({
            "type": "Feature",
            "properties": {"size_km": "5x5"},
            "geometry": mapping(t_wgs)
        })
        total_area += t.area / 1e6

    fc = {"type": "FeatureCollection", "features": feats}
    return GridResult(tiles_geojson=fc, n_tiles=len(feats), total_area_km2=total_area)

# ==========================
# Interface principal
# ==========================
st.title("📡 Planejador de Aquisições — GHGSat (5 km × 5 km)")
st.caption("Desenhe sua Área de Interesse, gere o fatiamento e exporte os footprints.")

col_left, col_mid, col_right = st.columns([1.15, 0.2, 1.15])

if "aoi_geojson" not in st.session_state:
    st.session_state["aoi_geojson"] = None
if "grid_result" not in st.session_state:
    st.session_state["grid_result"] = None

# ==========================
# Mapa da esquerda
# ==========================
with col_left:
    st.subheader("1) Desenhar área de interesse")
    m1 = folium.Map(location=[-15.78, -47.93], zoom_start=4, tiles="OpenStreetMap")

    from folium.plugins import Draw
    draw = Draw(export=False, filename='aoi.geojson', position='topleft',
                draw_options={
                    "polygon": True, "rectangle": True, "circle": False,
                    "polyline": False, "marker": False, "circlemarker": False
                },
                edit_options={"edit": True, "remove": True})
    draw.add_to(m1)

    left_map = st_folium(m1, height=550, returned_objects=["last_drawn_feature", "all_drawings"])

    if st.button("Confirmar AOI", type="primary"):
        feat = None
        if left_map and left_map.get("last_drawn_feature"):
            feat = left_map["last_drawn_feature"]
        elif left_map and _extract_features(left_map.get("all_drawings")):
            try:
                feats = _extract_features(left_map.get("all_drawings"))
                geoms = [shape(f["geometry"]) for f in feats if f.get("geometry")]
                if geoms:
                    merged = unary_union(geoms).buffer(0)
                    feat = {"type": "Feature", "geometry": mapping(merged), "properties": {}}
            except Exception as e:
                st.error(f"Erro ao ler desenhos: {e}")
        if not feat:
            st.warning("Desenhe um polígono antes de confirmar.")
        else:
            st.session_state["aoi_geojson"] = feat["geometry"]
            st.session_state["grid_result"] = None
            st.success("AOI confirmada com sucesso ✅")

# ==========================
# Botão planejar
# ==========================
with col_mid:
    st.write(" ")
    st.write(" ")
    plan_click = st.button("Planejar ➜", use_container_width=True)

# ==========================
# Mapa da direita
# ==========================
with col_right:
    st.subheader("2) Footprints gerados (5×5 km)")
    m2 = folium.Map(location=[-15.78, -47.93], zoom_start=4, tiles="OpenStreetMap")

    if plan_click:
        if not st.session_state["aoi_geojson"]:
            st.warning("Confirme uma AOI primeiro.")
        else:
            try:
                grid = generate_5km_grid_covering(st.session_state["aoi_geojson"])
                st.session_state["grid_result"] = grid
            except Exception as e:
                st.error(f"Erro no planejamento: {e}")

    if st.session_state["aoi_geojson"]:
        folium.GeoJson(
            {"type":"FeatureCollection","features":[{"type":"Feature","geometry":st.session_state["aoi_geojson"]}]},
            style_function=lambda x: {"fillColor":"#ffffff","color":"#000","weight":2,"fillOpacity":0.25}
        ).add_to(m2)

    if st.session_state["grid_result"]:
        folium.GeoJson(
            st.session_state["grid_result"].tiles_geojson,
            name="Footprints 5x5 km",
            style_function=lambda x: {"fillColor":"#cccccc","color":"#1f77b4","weight":2,"opacity":0.9,"fillOpacity":0.3}
        ).add_to(m2)
        st.metric("Número de footprints", st.session_state["grid_result"].n_tiles)
        st.metric("Área total (km²)", f'{st.session_state["grid_result"].total_area_km2:.2f}')
        geojson_str = json.dumps(st.session_state["grid_result"].tiles_geojson)
        st.download_button("Baixar GeoJSON", geojson_str, file_name="ghgsat_footprints.geojson", mime="application/geo+json")

    st_folium(m2, height=550)
