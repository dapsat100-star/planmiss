
# GHGSat 5×5 km Footprint Planner (Streamlit)

Planejador simples para gerar footprints quadrados de **5 km × 5 km** cobrindo uma **Área de Interesse (AOI)** desenhada pelo usuário.

## Como rodar localmente

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Como usar

1. Abra o app no navegador.
2. **Desenhe** um polígono (sua AOI) no mapa à esquerda e clique em **Confirmar AOI**.
3. Clique em **Planejar Aquisições**.
4. Veja à direita o **fatiamento em 5×5 km** e **baixe o GeoJSON** dos footprints.

## Notas técnicas

- A grade é gerada em **metros** num CRS **UTM** escolhido automaticamente pelo centróide da AOI (via `pyproj`).
- As geometrias são manipuladas com **Shapely** e os mapas com **Folium**.
- O resultado final é exibido em **WGS84 (EPSG:4326)**.
- Este repositório é um ponto de partida: adapte preços, políticas de priorização e integração com provedores.

## Deploy no Streamlit Cloud

1. Faça **fork** deste projeto no GitHub.
2. No Streamlit Cloud, aponte para o repositório e para o arquivo `app.py`.
3. Garanta que `requirements.txt` esteja presente.
