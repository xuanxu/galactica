import pytest
import galactica
from galactica import galaxy, settings


def test_default_initializing():
    g = galaxy.Galaxy()
    assert g.params == settings.default_settings()
