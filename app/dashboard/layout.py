"""
MST-Pipeline — Single-page layout (no sidebar navigation).
"""
from dash import Input, Output, dcc, html

from app.dashboard.app import app as dash_app


def create_layout():
    return html.Div(
        [
            dcc.Location(id="url", refresh=False),
            html.Div(id="page-content", className="p-4"),
        ]
    )


@dash_app.callback(Output("page-content", "children"), Input("url", "pathname"))
def render_page(pathname):
    from app.dashboard.pages.unified_page import get_layout
    return get_layout()
