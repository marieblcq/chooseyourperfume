import os
from streamlit.testing.v1 import AppTest

# Build the correct path to app_perfume.py
APP_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'app_perfume.py'))

def test_app_loads_without_errors():
    assert os.path.isfile(APP_PATH), f"Could not find: {APP_PATH}"
    app = AppTest.from_file(APP_PATH)
    app.run(timeout=15) 
    assert not app.exception, f"App raised unexpected exceptions: {app.exception}"


def test_app_title_displayed():
    app = AppTest.from_file(APP_PATH)
    app.run()
    title_elements = [e.value for e in app.markdown if "CHOOSE YOUR PERFUME" in e.value.upper()]
    assert title_elements, "Expected app title not found"


def test_categories_expanders_exist():
    app = AppTest.from_file(APP_PATH)
    app.run()
    categories = [e.label for e in app.expander]
    assert "Floral" in categories or "Woody" in categories or "Fruity" in categories, "Expected scent category expanders missing"


def test_surprise_me_button_exists():
    app = AppTest.from_file(APP_PATH)
    app.run()
    labels = [b.label for b in app.button]
    assert any("surprise" in label.lower() for label in labels), "Surprise button not found"
