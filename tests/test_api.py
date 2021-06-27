from psv_tools.api import func


def test_func():
    assert func('world') == 'Hello, world'
