"""
MST-Pipeline — Dash application initialization.
"""
import dash
import dash_bootstrap_components as dbc
import dash_uploader as du

from app.config import UPLOAD_DIR

app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.FLATLY],
    suppress_callback_exceptions=True,
    title="MST Pipeline",
    update_title="Loading...",
)

app.index_string = '''<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
        <style>body { background-color: #c4d7e8; }</style>
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>'''

server = app.server

# Allow large file uploads
server.config["MAX_CONTENT_LENGTH"] = 500 * 1024 * 1024  # 500 MB
server.config["MAX_FORM_MEMORY_SIZE"] = 500 * 1024 * 1024

# Configure dash-uploader for chunked FASTQ uploads
du.configure_upload(app, str(UPLOAD_DIR), use_upload_id=True)
