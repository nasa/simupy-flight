# Modifications

## Strict MathML comparison

Date: 2026-09-26

Originator: Sylvester Kaczmarek (github.com/sylvesterkaczmarek)

This Modification is derived from the Original Software, the NASA SimuPy Flight
Vehicle Toolkit. The originator consents to its characterization as a
Modification under the included NASA Open Source Agreement Version 1.3.

In `simupy_flight/parse_daveml.py`, translate the MathML `lt` operator using
SymPy's strict less-than relation, preserving false results at equality.
Add parser and generated-code regression tests for strict and non-strict
relations, scalar and array inputs, arithmetic operands, and invalid arity.
Other operators, numerical dynamics, constants, and dependencies are unchanged.
