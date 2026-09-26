"""Exercise relational MathML through parsing and generated NumPy code."""

import importlib.util
import operator

from lxml import etree
import numpy as np
import pytest
import sympy as sp

from simupy_flight.parse_daveml import ProcessDaveML


OPERATIONS = {"lt": operator.lt, "leq": operator.le,
              "gt": operator.gt, "geq": operator.ge}
SYMBOLIC = {"lt": sp.Lt, "leq": sp.Le, "gt": sp.Gt, "geq": sp.Ge}
MATHML = "http://www.w3.org/1998/Math/MathML"


def _generate(tmp_path, expression):
    source = tmp_path / "relation.dml"
    source.write_text(
        '<DAVEfunc xmlns="http://daveml.org/2010/DAVEML" '
        'xmlns:m="' + MATHML + '">'
        '<variableDef varID="x"><isInput/></variableDef>'
        '<variableDef varID="y"><isInput/></variableDef>'
        '<variableDef varID="result"><isOutput/><calculation><m:math>'
        + expression + '</m:math></calculation></variableDef></DAVEfunc>',
        encoding="utf-8",
    )
    parser = ProcessDaveML(str(source), outdir=str(tmp_path))
    spec = importlib.util.spec_from_file_location(
        "generated_relation", tmp_path / "relation.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return parser, module.relation


@pytest.fixture(params=OPERATIONS)
def relation(request, tmp_path):
    tag = request.param
    parser, function = _generate(
        tmp_path, '<m:apply><m:' + tag
        + '/><m:ci>x</m:ci><m:ci>y</m:ci></m:apply>'
    )
    return tag, parser, function


@pytest.mark.parametrize("left,right", [
    (-2.0, -2.0), (0.0, 0.0), (2.0, 2.0),
    (-1.0, 1.0), (1.0, -1.0), (1.0, np.nextafter(1.0, 2.0)),
])
def test_generated_relation_matches_operator(relation, left, right):
    tag, _, function = relation
    actual = function(left, right)
    assert actual.shape == (1,)
    assert bool(actual[0]) is bool(OPERATIONS[tag](left, right))


def test_symbolic_relation_retains_strictness(relation):
    tag, parser, _ = relation
    actual = parser.extra_assignments[parser.symbol_table["result"]]
    expected = SYMBOLIC[tag](
        parser.symbol_table["x"], parser.symbol_table["y"]
    )
    assert actual == expected


def test_generated_relation_accepts_array_inputs(relation):
    tag, _, function = relation
    left = np.array([-2.0, 0.0, 2.0, 1.0, 1.0])
    right = np.array([-2.0, 0.0, 2.0, 2.0, 0.0])
    np.testing.assert_array_equal(function(left, right)[0],
                                  OPERATIONS[tag](left, right))


@pytest.mark.parametrize("count", [0, 1, 3])
def test_relations_still_require_two_operands(relation, count):
    tag, parser, _ = relation
    node = etree.fromstring(
        '<apply xmlns="' + MATHML + '"><' + tag + '/>'
        + '<ci>x</ci>' * count + '</apply>'
    )
    with pytest.raises(ValueError, match="binary operation"):
        parser.xml_to_sympy_expr(node)


def test_strict_relation_with_arithmetic_operand(tmp_path):
    _, function = _generate(
        tmp_path,
        '<m:apply><m:lt/><m:apply><m:plus/><m:ci>x</m:ci><m:cn>1</m:cn>'
        '</m:apply><m:ci>y</m:ci></m:apply>',
    )
    assert not function(2.0, 3.0)[0]
    assert function(2.0, 4.0)[0]
