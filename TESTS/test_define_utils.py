from define_utils import add, is_even, capitalize_words

def test_add():
    assert add(2, 3) == 5
    assert add(-1, 1) == 0
    assert add(0, 0) == 0
    pass

def test_is_even():
    assert is_even(2) is True
    assert is_even(3) is False
    assert is_even(0) is True
    pass


def test_capitalize_words():
    assert capitalize_words("hello world") == "Hello World"
    assert capitalize_words("python is fun") == "Python Is Fun"
    assert capitalize_words("") == ""
    pass
