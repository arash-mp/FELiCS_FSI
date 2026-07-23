#  ___________________________________   _______________________________________________________
# /-----------------------------------\ /-------------------------------------------------------\
# |   (         (                (     |  This source code is part of FELiCS                     |
# |   )\ )     ) )        (      )\ )  |  (F)inite (E)lement (Li)nearized (C)ombustion (S)olver  |
# |  (()/(  (  (()/( (    )\   (()/(   |                                                         |
# |  /(_)) )\  /(_)))\  (((_)  /(_))   |  Licensed under the GNU GPLv3                           |
# |  (_)_)((_) (_)) ((_) )\___ (_))    |                                                         |
# |  | __|| __|| |   (_)((/ __|/ __|   |  (C) 2018-2025: The FELiCS Developers (www.felics.eu)   |
# |  | _| | _| | |__ | | | (__ \__ \   |  Visit          www.felics.eu                           |
# |  |_|  |___||____||_|  \___||___/   |  Contact        info@felics.eu                          |
# \___________________________________/ \_______________________________________________________/
#

import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import logging

from FELiCS.Misc.logging import CustomFormatter


def test_custom_formatter_includes_record_fields_in_non_debug_mode():
    formatter = CustomFormatter(debug_mode=False)
    record = logging.LogRecord(
        name="felics",
        level=logging.INFO,
        pathname=__file__,
        lineno=42,
        msg="hello world",
        args=None,
        exc_info=None,
        func="my_func",
    )

    formatted = formatter.format(record)

    assert "hello world" in formatted
    assert "my_func" in formatted
    assert "42" in formatted


def test_custom_formatter_includes_record_fields_in_debug_mode():
    formatter = CustomFormatter(debug_mode=True)
    record = logging.LogRecord(
        name="felics",
        level=logging.DEBUG,
        pathname=__file__,
        lineno=42,
        msg="hello world",
        args=None,
        exc_info=None,
        func="my_func",
    )

    formatted = formatter.format(record)

    assert "hello world" in formatted
    assert "my_func" in formatted
    assert "42" in formatted
