import alpaqa as pa

def test_version():
    assert pa.__version__ == pa.__c_version__
