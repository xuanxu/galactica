import pytest
import galactica


def test_version():
    assert galactica.__version__ is not None
