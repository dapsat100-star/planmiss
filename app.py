# -*- coding: utf-8 -*-
# app.py — GHGSat 5x5 km Footprint Planner (Streamlit + Folium)
#
# Como usar (resumo):
# 1) Desenhe um polígono (Área de Interesse) no mapa à esquerda e clique em "Confirmar AOI".
# 2) Clique em "Planejar Aquisições" para gerar os footprints 5 km x 5 km que cobrem a AOI.
# 3) Veja o fatiamento à direita e, se quiser, exporte os footprints em GeoJSON.
#
# Observações:
# - O cálculo dos quadrados de 5 km é feito em metros em um sistema projetado local (UTM),
#   derivado automaticamente do centróide da AOI. Depois, convertemos de volta para WGS84.
# - Dependências: streamlit, folium, streamlit-folium, shapely, pyproj

import json
from typing import List, Tuple
from dataclasses import dataclass

import streamlit as st
import folium
from streamlit_folium import st_folium
from shapely.geometry import Polygon, shape, box, mapping
from shapely.ops import unary_union
from shapely import wkt
from pyproj import CRS, Transformer


        # Tenta extrair todas as feições desenhadas de forma tolerante
def _extract_features(all_drawings_obj):
            feats = []
            if not all_drawings_obj:
                return feats
            # Caso 1: dicionário estilo GeoJSON FeatureCollection
            if isinstance(all_drawings_obj, dict):
                if all_drawings_obj.get("type") == "FeatureCollection" and all_drawings_obj.get("features"):
                    return list(all_drawings_obj["features"])
                # alguns retornam {"features":[...]} sem o "type"
                if all_drawings_obj.get("features"):
                    return list(all_drawings_obj["features"])
            # Caso 2: lista de feições
            if isinstance(all_drawings_obj, list):
                for item in all_drawings_obj:
                    if isinstance(item, dict) and item.get("geometry"):
                        feats.append(item)
            return feats


st.set_page_config(page_title="Planejador GHGSat 5x5 km", layout="wide")

# -------------------------- Utils de Projeção --------------------------

def utm_crs_from_lonlat(lon: float, lat: float) -> CRS:
    zone = int((lon + 180) // 6) + 1
    south = lat < 0
    return CRS.from_dict({
        "proj": "utm",
        "zone": zone,
        "south": south
    })

def to_proj(geom_geojson, crs_to: CRS):
    """Converte GeoJSON (WGS84) para shapely em CRS projetado."""
    geom = shape(geom_geojson)
    transformer = Transformer.from_crs("epsg:4326", crs_to, always_xy=True)
    def _reproj_coords(coords):
        x, y = transformer.transform(coords[0], coords[1])
        return (x, y)
    if geom.geom_type == "Polygon":
        exterior = [ _reproj_coords(c) for c in list(geom.exterior.coords) ]
        interiors = []
        for ring in geom.interiors:
            interiors.append([ _reproj_coords(c) for c in list(ring.coords) ])
        return Polygon(exterior, interiors)
    elif geom.geom_type == "MultiPolygon":
        # Unifica em um polígono (buffer(0) para limpar geometrias)
        polys = []
        for g in geom.geoms:
            exterior = [ _reproj_coords(c) for c in list(g.exterior.coords) ]
            interiors = []
            for ring in g.interiors:
                interiors.append([ _reproj_coords(c) for c in list(ring.coords) ])
            polys.append(Polygon(exterior, interiors))
        return unary_union(polys).buffer(0)
    else:
        raise ValueError("Somente Polygon/MultiPolygon suportados para AOI.")

def to_wgs84(geom: Polygon, crs_from: CRS):
    transformer = Transformer.from_crs(crs_from, "epsg:4326", always_xy=True)
    def _reproj_coords(coords):
        x, y = transformer.transform(coords[0], coords[1])
        return (x, y)
    exterior = [ _reproj_coords(c) for c in list(geom.exterior.coords) ]
    interiors = []
    for ring in geom.interiors:
        interiors.append([ _reproj_coords(c) for c in list(ring.coords) ])
    return Polygon(exterior, interiors)

# -------------------------- Grid 5x5 km --------------------------

@dataclass
class GridResult:
    tiles_geojson: dict
    n_tiles: int
    total_area_km2: float

def generate_5km_grid_covering(aoi_geojson: dict) -> GridResult:
    # Projeta AOI para um CRS em metros (UTM)
    centroid = shape(aoi_geojson).centroid
    crs_utm = utm_crs_from_lonlat(centroid.x, centroid.y)
    aoi_m = to_proj(aoi_geojson, crs_utm)

    # Cria bbox e grade 5x5 km cobrindo a AOI
    cell = 5000.0  # metros
    minx, miny, maxx, maxy = aoi_m.bounds
    # margem de 1 célula
    minx -= cell; miny -= cell; maxx += cell; maxy += cell

    # Alinha à malha 5 km (para reprodutibilidade)
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

    # Converte de volta para WGS84 e constrói GeoJSON FeatureCollection
    feats = []
    total_area = 0.0
    for t in tiles:
        if t.is_empty:
            continue
        # área em CRS projetado (m²)
        total_area += float(t.area) / 1e6  # km²
        # converte para WGS84 para visualização/exportação
        t_wgs = to_wgs84(t, crs_utm)
        feats.append({
            "type": "Feature",
            "properties": {"size_km": "5x5"},
            "geometry": mapping(t_wgs)
        })

    fc = {"type": "FeatureCollection", "features": feats}
    return GridResult(tiles_geojson=fc, n_tiles=len(feats), total_area_km2=total_area)

# -------------------------- UI --------------------------

st.title("Planejador de Aquisições — GHGSat (5 km x 5 km)")
st.caption("Desenhe sua Área de Interesse (AOI), gere o fatiamento em footprints 5×5 km e exporte o plano.")

col_left, col_mid, col_right = st.columns([1.15, 0.2, 1.15])

# Sessão
if "aoi_geojson" not in st.session_state:
    st.session_state["aoi_geojson"] = None
if "grid_result" not in st.session_state:
    st.session_state["grid_result"] = None

with col_left:
    st.subheader("1) Desenhar área de interesse")
    st.write("Use a ferramenta de desenho no mapa. Ao terminar, clique em **Confirmar AOI**.")
    # Mapa base
    m1 = folium.Map(location=[-15.78, -47.93], zoom_start=4, tiles="OpenStreetMap")
    # Plug-in de desenho
    from folium.plugins import Draw
    draw = Draw(export=False, filename='aoi.geojson', position='topleft',
                draw_options={
                    "polygon": True, "rectangle": True, "circle": False,
                    "polyline": False, "marker": False, "circlemarker": False
                },
                edit_options={"edit": True, "remove": True})
    draw.add_to(m1)

    left_map = st_folium(m1, height=550, width=None, returned_objects=["last_drawn_feature", "all_drawings"])

    # Botão confirmar
    if st.button("Confirmar AOI", type="primary"):
        # Prioriza última feição desenhada; se não, tenta pegar o conjunto completo
        feat = None

        # 1) Última feição desenhada (quando disponível)
        if isinstance(left_map, dict) and left_map.get("last_drawn_feature"):
            feat = left_map["last_drawn_feature"]

        # 2) Caso contrário, tenta consolidar "all_drawings" com múltiplas feições
        if feat is None and isinstance(left_map, dict):
            drawings = left_map.get("all_drawings")
            features = None
            if isinstance(drawings, dict) and "features" in drawings:
                features = drawings.get("features")
            elif isinstance(drawings, list):
                # algumas versões retornam lista de features diretamente
                features = drawings
            if features:
                try:
                    geoms = []
                    for f in features:
                        geom = f.get("geometry") if isinstance(f, dict) else None
                        if not geom:
                            continue
                        gtype = geom.get("type")
                        if gtype in ("Polygon", "MultiPolygon"):
                            geoms.append(shape(geom))
                        elif gtype == "GeometryCollection" and geom.get("geometries"):
                            for g in geom["geometries"]:
                                if g.get("type") in ("Polygon", "MultiPolygon"):
                                    geoms.append(shape(g))
                    if geoms:
                        merged = unary_union(geoms).buffer(0)
                        feat = {"type":"Feature","properties":{},"geometry":mapping(merged)}
                except Exception as e:
                    st.error(f"Não consegui ler os desenhos ({e}). Tente novamente com um único polígono.")

        if not feat:
            st.warning("Desenhe um polígono antes de confirmar.")
        else:
            st.session_state["aoi_geojson"] = feat["geometry"]
            st.session_state["grid_result"] = None
            st.success("AOI confirmada! Agora clique em **Planejar Aquisições**.")

with col_mid:
    st.write(" ")
    st.write(" ")
    st.write(" ")
    plan_click = st.button("Planejar Aquisições ➜", use_container_width=True)

with col_right:
    st.subheader("2) Fatiamento do Satélite (5×5 km)")
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

    # Desenha AOI (se existir)
    if st.session_state["aoi_geojson"]:
        folium.GeoJson(
            {"type":"FeatureCollection","features":[{"type":"Feature","properties":{"name":"AOI"},"geometry":st.session_state["aoi_geojson"]}]},
            name="AOI",
            style_function=lambda x: {"fillColor":"#ffffff","color":"#000000","weight":2,"opacity":0.8,"fillOpacity":0.25}
        ).add_to(m2)

    # Desenha grid (se existir)
    if st.session_state["grid_result"]:
        folium.GeoJson(
            st.session_state["grid_result"].tiles_geojson,
            name="Footprints 5x5 km",
            style_function=lambda x: {"fillColor":"#cccccc","color":"#1f77b4","weight":2,"opacity":0.9,"fillOpacity":0.3},
            tooltip=folium.GeoJsonTooltip(fields=[], aliases=[], sticky=False)
        ).add_to(m2)
        st.metric("Número de footprints (5×5 km)", st.session_state["grid_result"].n_tiles)
        st.metric("Área total estimada (km²)", f'{st.session_state["grid_result"].total_area_km2:.2f}')

        # Exportar
        geojson_str = json.dumps(st.session_state["grid_result"].tiles_geojson, ensure_ascii=False)
        st.download_button("Baixar Footprints (GeoJSON)", data=geojson_str, file_name="ghgsat_footprints_5x5.geojson", mime="application/geo+json")

    st_folium(m2, height=550, width=None)

st.markdown("---")
with st.expander("Ajuda e Observações"):
    st.markdown("""
- **Dica**: Desenhe um *retângulo* se quiser um plano mais compacto; o algoritmo sempre cobre toda a AOI.
- **Precisão**: As áreas e os limites são aproximados (dependem da projeção local escolhida automaticamente).
- **Custo**: Você pode estender este app para estimar custo por footprint (ex.: valor por cena), somando `n_tiles × preço`.
- **Integrações**: Para exportar para QGIS, use o GeoJSON baixado.
    """)
 
