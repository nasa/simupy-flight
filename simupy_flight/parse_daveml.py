import os
import re
from lxml import etree
import sympy as sp, numpy as np, ndsplines
import warnings
from simupy import codegen

delimeter_regex = re.compile("\s+|(?<!\d)[.](?!\d)|,")


def isTag(node, tagsuffix=""):
    return (
        not isinstance(node, etree._Comment)
        and isinstance(node.tag, str)
        and node.tag.endswith(tagsuffix)
    )


def filterTags(nodes, tag=""):
    return [node for node in nodes if isTag(node, tag)]


Clip = sp.Function("clip")

run_check_code_template = """
def run_checks():
    for check_input, check_output in check_data:
        if not numpy.allclose( {funcname}(*check_input), check_output):
            raise ValueError("Check for {funcname} failed!")
    print("All checks for {funcname} passed.")

if __name__ == "__main__":
    run_checks()
"""


class ProcessDaveML:
    def __init__(self, filename, outdir="."):
        self.input_vars = []
        self.extra_assignments = {}
        self.output_vars = []
        self.symbol_table = {}
        self.breakpoints = {}
        self.constants = {}
        self.splines = {}
        self.extra_tables = {}
        self.check_data = []

        with open(filename, "r") as xmlfile:
            text_data = xmlfile.read().replace("\n", "")

        xml_parsed_data = etree.XML(text_data)
        self.process_xml(xml_parsed_data)

        path, fname = os.path.split(filename)
        funcname = fname.split(".")[0]
        # mkdir? and add __init__.py
        # no, it could just be one file with everything
        # will need to dump spline data, so could do that as a dir or just as names
        # also need to save constants
        constant_modules = {constant: "" for constant in self.constants}
        function_modules = {spline.func: "" for spline in self.splines}
        MNPP = codegen.ModuleNumPyPrinter(
            constant_modules=constant_modules,
            function_modules=function_modules,
        )

        pre_lines = []
        if self.splines:
            pre_lines.append("import ndsplines")
        pre_lines.append("")

        if self.constants:
            pre_lines.append("")
            for var, val in self.constants.items():
                pre_lines.append("%s = %e" % (str(var), val))

        array_printer = codegen.ArrayNumPyPrinter()

        if self.splines:
            pre_lines.append("\n")
            for var, spline in self.splines.items():
                base_name = str(var.func)
                knot_vars = []
                for knot_idx, knot_arr in enumerate(spline.knots):
                    knot_var = "%s_knots_%d" % (base_name, knot_idx)
                    knot_value = array_printer.doprint(knot_arr)
                    knot_vars.append(knot_var)
                    pre_lines.append("%s = %s" % (knot_var, knot_value))

                knots_to_init_str = "[%s]" % ",".join(
                    [knot_var for knot_var in knot_vars]
                )

                to_save_coeff = spline.coefficients.reshape(
                    spline.coefficients.shape[: spline.xdim] + spline.yshape
                )
                new_line_str = "]" * (to_save_coeff.ndim - 1) + ","
                coefficient = array_printer.doprint(to_save_coeff).replace(
                    new_line_str, new_line_str + "\n"
                )
                coeff_varname = "%s_coeffs" % (base_name)
                pre_lines.append("%s = %s" % (coeff_varname, coefficient))
                pre_lines.append(
                    "%s = ndsplines.NDSpline(%s, %s, 1)\n"
                    % (base_name, knots_to_init_str, coeff_varname)
                )

        pre_line_str = "\n".join(pre_lines)
        printer = codegen.ModulePrinter(MNPP)
        if self.input_vars:
            code = printer.codeprint(
                [
                    [
                        funcname,
                        sp.Array(self.input_vars),
                        sp.Array(self.output_vars),
                        self.extra_assignments,
                    ]
                ],
                extra_pre_lines=pre_line_str,
            )
        else:
            code = "import numpy\n" + pre_line_str

        if self.check_data:
            check_data_code = "\ncheck_data = [\n [%s] \n]\n\n" % "],\n\n [".join(
                [
                    "%s,\n %s"
                    % (list(check_input.values()), list(check_output.values()))
                    for check_input, check_output in self.check_data
                ]
            )
            run_check_code = run_check_code_template.format(funcname=funcname)
            code = code + check_data_code + run_check_code

        out_fp = os.path.join(outdir, f"{funcname}.py")
        with open(out_fp, "w") as codefile:
            codefile.write(code)

    def xml_to_sympy_expr(self, xml_node):
        if isTag(xml_node, "/MathML}cn"):
            return float(xml_node.text)
        elif isTag(xml_node, "/MathML}ci"):
            if xml_node.text not in self.symbol_table:
                self.symbol_table[xml_node.text] = sp.Symbol(xml_node.text)
                self.extra_assignments[self.symbol_table[xml_node.text]] = None
            return self.symbol_table[xml_node.text]

        elif not isTag(xml_node, "/MathML}apply"):
            raise ValueError("Do not know hot to parse node %s" % xml_node)

        try:
            func, *arg_children = filterTags(xml_node.getchildren())
        except ValueError:
            raise
        processed_args = []
        for child in arg_children:
            processed_args.append(self.xml_to_sympy_expr(child))

        if isTag(func, "/MathML}plus"):
            out_expr = sp.Add(*processed_args)
        elif isTag(func, "/MathML}times"):
            out_expr = sp.Mul(*processed_args)
        elif isTag(func, "/MathML}power"):
            if len(processed_args) != 2:
                raise ValueError("Power is a binary operation!")
            out_expr = sp.Pow(*processed_args)
        elif isTag(func, "/MathML}divide"):
            if len(processed_args) != 2:
                raise ValueError("Division is a binary operation!")
            out_expr = sp.Mul(processed_args[0], sp.Pow(processed_args[1], -1))
        elif isTag(func, "/MathML}minus"):
            if (len(processed_args) > 2) or (len(processed_args) == 0):
                raise ValueError("Subtraction is a unitary or binary operation!")
            elif len(processed_args) == 2:
                out_expr = sp.Add(processed_args[0], sp.Mul(processed_args[1], -1))
            elif len(processed_args) == 1:
                out_expr = sp.Mul(processed_args[0], -1)
        elif isTag(func, "/MathML}abs"):
            if len(processed_args) != 1:
                raise ValueError("Absolute value is a unitary operation!")
            out_expr = sp.Abs(processed_args[0])

        elif isTag(func, "/MathML}cos"):
            if len(processed_args) != 1:
                raise ValueError("cos is a unitary operation!")
            out_expr = sp.cos(processed_args[0])
        elif isTag(func, "/MathML}sin"):
            if len(processed_args) != 1:
                raise ValueError("sin is a unitary operation!")
            out_expr = sp.sin(processed_args[0])
        elif isTag(func, "/MathML}csymbol"):

            if func.text == "atan2":
                if len(processed_args) != 2:
                    raise ValueError("atan2 is a binary operation!")
                out_expr = sp.atan2(*processed_args)
            else:
                raise ValueError("Unknown csymbol %s" % func.text)

        elif isTag(func, "/MathML}piecewise"):
            *pieces, otherwise = func.getchildren()
            otherwise_children = filterTags(otherwise.getchildren())
            conditions = []
            values = []
            if len(otherwise_children) != 1:
                raise ValueError("Too many `otherwise` children")
            for piece in filterTags(pieces):
                value_child, condition_child = filterTags(piece.getchildren())
                values.append(self.xml_to_sympy_expr(value_child))
                conditions.append(self.xml_to_sympy_expr(condition_child))

            values.append(self.xml_to_sympy_expr(otherwise_children[0]))
            conditions.append(True)
            out_expr = sp.Piecewise(
                *[(value, condition) for value, condition in zip(values, conditions)]
            )
        elif isTag(func, "/MathML}lt"):
            if len(processed_args) != 2:
                raise ValueError("Less than is a binary operation!")
            out_expr = sp.Le(*processed_args)
        elif isTag(func, "/MathML}gt"):
            if len(processed_args) != 2:
                raise ValueError("Greater than is a binary operation!")
            out_expr = sp.Gt(*processed_args)
        elif isTag(func, "/MathML}leq"):
            if len(processed_args) != 2:
                raise ValueError("Less than equal to is a binary operation!")
            out_expr = sp.Le(*processed_args)
        elif isTag(func, "/MathML}geq"):
            if len(processed_args) != 2:
                raise ValueError("Greater than equal to is a binary operation!")
            out_expr = sp.Ge(*processed_args)

        else:
            raise ValueError("Unknown func tag, %s" % func.tag)

        if out_expr is None:
            raise ValueError("Calculated None")
        else:
            return out_expr

    def check_dependencies(self, atoms):
        for symbol in atoms:
            if (
                (symbol not in self.input_vars)
                and (symbol not in self.extra_assignments.keys())
                and (symbol not in self.constants)
            ):
                raise ValueError("Unspecified dependency")

    def process_xml(self, xml_parsed_data):
        for idx, child in enumerate(xml_parsed_data.getchildren()):
            if isinstance(child.tag, str):
                if isTag(child, "/DAVEML}variableDef"):
                    if child.attrib["varID"] in self.symbol_table:
                        warnings.warn(
                            "%s variable already in symbol table!"
                            % child.attrib["varID"]
                        )
                    self.symbol_table[child.attrib["varID"]] = sp.Symbol(
                        child.attrib["varID"]
                    )

                    hasLimits = ("minValue" in child.attrib) or (
                        "maxValue" in child.attrib
                    )

                    isInput = False
                    isOutput = False
                    isAssigned = False

                    for varChild in child.getchildren():
                        if isTag(varChild, "/DAVEML}isInput"):
                            isInput = True
                        if isTag(varChild, "/DAVEML}isOutput"):
                            isOutput = True
                        if isTag(varChild, "/DAVEML}calculation"):
                            mathTags = filterTags(varChild.getchildren())
                            if len(mathTags) != 1:
                                raise ValueError(
                                    "A calculation tag should hold only one non-comment tag"
                                )
                            if not isTag(mathTags[0], "/MathML}math"):
                                raise ValueError(
                                    "A calculation tag should hold only one math tag"
                                )

                            for mathChild in mathTags[
                                0
                            ].getchildren():  # iterate over children of math tag
                                processed_apply = False
                                if isinstance(mathChild, etree._Comment):
                                    pass
                                elif isTag(mathChild, "/MathML}apply"):
                                    if processed_apply:
                                        raise ValueError(
                                            "Trying to apply multiple times within a Math tag"
                                        )
                                    processed_apply = True
                                    expr = self.xml_to_sympy_expr(mathChild)
                                    self.check_dependencies(expr.atoms(sp.Symbol))
                                    self.extra_assignments[
                                        self.symbol_table[child.attrib["varID"]]
                                    ] = expr
                                    isAssigned = True
                                else:
                                    raise ValueError(
                                        "Unexpected tag under calculation.math %s"
                                        % mathChild
                                    )
                    if ("initialValue" in child.attrib.keys()) and not isInput:
                        # we have no way of knowing if this is a constant or an internal variable with "initial value" -- not sure why this needs that!
                        self.constants[
                            self.symbol_table[child.attrib["varID"]]
                        ] = float(child.attrib["initialValue"])
                    if isInput:
                        self.input_vars.append(self.symbol_table[child.attrib["varID"]])
                    if isOutput:
                        self.output_vars.append(
                            self.symbol_table[child.attrib["varID"]]
                        )
                    if not isAssigned:
                        # declared with default value (may be overwritten) OR unspecified
                        self.extra_assignments[
                            self.symbol_table[child.attrib["varID"]]
                        ] = None
                    elif hasLimits:  # isAssigned
                        self.extra_assignments[
                            self.symbol_table[child.attrib["varID"]]
                        ] = sp.Function("numpy.clip")(
                            self.extra_assignments[
                                self.symbol_table[child.attrib["varID"]]
                            ],
                            float(child.attrib.get("minValue", -np.inf)),
                            float(child.attrib.get("maxValue", np.inf)),
                        )

                elif isTag(child, "/DAVEML}fileHeader"):
                    print("Skipping header")

                elif isTag(child, "/DAVEML}breakpointDef"):
                    description, bpvals = child.getchildren()
                    self.breakpoints[child.attrib["bpID"]] = np.array(
                        [float(num) for num in bpvals.text.split(",")]
                    )

                elif isTag(child, "/DAVEML}function") or isTag(
                    child, "/DAVEML}griddedTableDef"
                ):
                    # function
                    #   independentVarRef
                    #   dependentVarRef
                    #   functionDefn
                    #     griddedTableDef gtID
                    #       breakpoint_refs
                    #       datatable
                    input_vars = []
                    # DaveML assumes only 1 output for function

                    from_function = isTag(child, "/DAVEML}function")
                    from_griddedTable = isTag(child, "/DAVEML}griddedTableDef")
                    make_spline = True
                    for functionChild in child.getchildren():
                        if isinstance(functionChild, etree._Comment):
                            pass
                        elif from_function and isTag(
                            functionChild, "/DAVEML}independentVarRef"
                        ):
                            input_vars.append(
                                self.symbol_table[functionChild.attrib["varID"]]
                            )
                        elif from_function and isTag(
                            functionChild, "/DAVEML}dependentVarRef"
                        ):
                            output_var = self.symbol_table[
                                functionChild.attrib["varID"]
                            ]
                        elif (
                            from_function
                            and isTag(functionChild, "/DAVEML}functionDefn")
                        ) or from_griddedTable:
                            if from_function:
                                functionDefnChildren = filterTags(
                                    functionChild.getchildren()
                                )
                                if len(functionDefnChildren) != 1:
                                    raise ValueError(
                                        "A functionDefn tag should hold only one non-comment tag"
                                    )

                                if isTag(
                                    functionDefnChildren[0], "/DAVEML}griddedTableDef"
                                ):
                                    griddedTableDef = functionDefnChildren[0]
                                elif isTag(
                                    functionDefnChildren[0], "/DAVEML}griddedTableRef"
                                ):
                                    spline = self.extra_tables[
                                        functionDefnChildren[0].attrib["gtID"]
                                    ]
                                    make_spline = False

                                else:
                                    raise ValueError(
                                        "A functionDefn tag should hold only one griddedTableDef or griddedTableRef tag"
                                    )
                            elif from_griddedTable:
                                griddedTableDef = child

                            if make_spline:
                                breakpoint_refs, data_table = None, None
                                for griddedTableDefChild in filterTags(
                                    griddedTableDef.getchildren()
                                ):
                                    if isTag(
                                        griddedTableDefChild, "/DAVEML}breakpointRefs"
                                    ):
                                        if breakpoint_refs is None:
                                            breakpoint_refs = griddedTableDefChild
                                        else:
                                            raise ValueError(
                                                "Too many breakpoint_refs tags"
                                            )

                                    if isTag(griddedTableDefChild, "/DAVEML}dataTable"):
                                        if data_table is None:
                                            data_table = griddedTableDefChild
                                        else:
                                            raise ValueError("Too many data_table tags")

                                knots_arrays = []
                                for breakpoint_ref in breakpoint_refs:
                                    knots_arrays.append(
                                        self.breakpoints[breakpoint_ref.attrib["bpID"]]
                                    )

                                data_table = np.array(
                                    [
                                        float(num)
                                        for num in delimeter_regex.split(
                                            data_table.xpath("string()").strip()
                                        )
                                        if num != ""
                                    ]
                                ).reshape(
                                    [knots_array.size for knots_array in knots_arrays]
                                )
                                spline = ndsplines.make_interp_spline(
                                    knots_arrays, data_table, degrees=1
                                )

                                if from_griddedTable:
                                    self.extra_tables[
                                        griddedTableDef.attrib["gtID"]
                                    ] = spline
                                    break

                            if from_function:
                                spline_sym_func = sp.symbols(
                                    functionChild.attrib["name"], cls=sp.Function
                                )(tuple(input_vars))
                                spline_sym_func.shape = (1,)
                                self.splines[spline_sym_func] = spline
                                self.check_dependencies(input_vars)
                                # self.extra_assignments[output_var] = spline_sym_func
                                self.extra_assignments[output_var] = spline_sym_func

                        elif isTag(functionChild, "/DAVEML}description") or isTag(
                            functionChild, "/DAVEML}provenance"
                        ):
                            pass
                        else:
                            raise ValueError(
                                "Unexpected child %s of function %s "
                                % (functionChild.tag, child.tag)
                            )

                elif isTag(child, "/DAVEML}checkData"):
                    for checkDataTag in filterTags(child.getchildren()):
                        inputs = {}
                        outputs = {}

                        # TODO: defensive checks for tags? -- this is a common pattern, I could definitely refactor
                        # maybe lxml/ElementTree has a way to do it?
                        inputTag = None
                        outputTag = None
                        for checkTag in filterTags(checkDataTag.getchildren()):
                            if isTag(checkTag, "/DAVEML}checkInputs"):
                                if inputTag is not None:
                                    raise ValueError("Too many checkInputs tags")
                                inputTag = checkTag
                            elif isTag(checkTag, "/DAVEML}checkOutputs"):
                                if outputTag is not None:
                                    raise ValueError("Too many checkOutputs tags")
                                outputTag = checkTag
                            elif isTag(checkTag, "/DAVEML}internalValues"):
                                pass
                            else:
                                raise ValueError("Unsupported tag")

                        for local_dict, global_dict, tag in zip(
                            [inputs, outputs],
                            [self.input_vars, self.output_vars],
                            [inputTag, outputTag],
                        ):
                            for signal in tag:
                                # TODO: unit conversion check? Or document assumption of consistent units
                                # to follow the spec, would need more careful handling... https://daveml.org/DTDs/1p7b1/Ref/re40.html
                                signalName = None
                                signalUnits = None
                                signalValue = None
                                signalTol = None

                                for signalChild in filterTags(signal.getchildren()):
                                    if isTag(signalChild, "/DAVEML}signalName"):
                                        if signalName is not None:
                                            raise ValueError("Too many signalName tags")
                                        signalName = signalChild
                                    elif isTag(signalChild, "/DAVEML}signalID"):
                                        raise ValueError("signalID tag unsupported")
                                    elif isTag(signalChild, "/DAVEML}signalUnits"):
                                        if signalUnits is not None:
                                            raise ValueError(
                                                "Too many signalUnits tags"
                                            )
                                        signalUnits = signalChild
                                    elif isTag(signalChild, "/DAVEML}signalValue"):
                                        if signalValue is not None:
                                            raise ValueError(
                                                "Too many signalValue tags"
                                            )
                                        signalValue = signalChild
                                    elif isTag(signalChild, "/DAVEML}tol"):
                                        pass
                                    else:
                                        raise ValueError("Unsupported tag")
                                # TODO: input/output use "name" here instead of "varID", so to validate the data would require tracking it
                                # aero example preserves order, so we will assume that as well...
                                # TODO: document assumption

                                # if self.symbol_table[signalName.text] not in global_dict:
                                #     raise ValueError("Unrecognized variable %s" % signalName.text)
                                local_dict[signalName.text] = float(signalValue.text)
                        self.check_data.append([inputs, outputs])

                else:
                    raise ValueError("Unknown tag: %s" % child.tag)
            else:
                if not isinstance(child, etree._Comment):
                    print("not str:", idx, child.tag, child.attrib, child.text)

        pruned_extra_assignments = {}
        for assigned_var in self.extra_assignments.keys():
            if self.extra_assignments[assigned_var] is not None:
                # if it was assigned keep it
                if hasattr(self.extra_assignments[assigned_var], "shape"):
                    if self.extra_assignments[assigned_var].shape == (1,):
                        pruned_extra_assignments[
                            sp.Array([assigned_var])
                        ] = self.extra_assignments[assigned_var]
                    else:
                        # TODO: need to figure out how to define a single symbol as a matrix/array with shape
                        # TODO: and whatever that syntax is, it needs to make sure this works/has a fix
                        raise ValueError(
                            "Cannot handle assignment of non-scalar shaped functions'"
                        )
                else:
                    pruned_extra_assignments[assigned_var] = self.extra_assignments[
                        assigned_var
                    ]
                # and remove from constants if present
                self.constants.pop(assigned_var, None)
            elif (
                assigned_var not in self.constants
                and assigned_var not in self.input_vars
            ):  # not assigned and not in constant
                warnings.warn("%s was never assigned, treating as input" % assigned_var)
                self.input_vars.append(assigned_var)
            # if unassigned but in constants, don't add to pruned so it is fetched from constants
        self.extra_assignments = pruned_extra_assignments
        print("Completed")
