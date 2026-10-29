"""Static audit of top-level scalar readers; no production request introspection.

Nested objects/arrays and library ``from_dict`` schemas are outside this
inventory. Variable-key comprehensions are bound explicitly to their source
default maps; unresolved reader sites are reported rather than discarded.
"""

import ast
import inspect
import textwrap

from rocket_propulsion.api import transient_routes as transient

NUMBER_READERS = frozenset({'_number', '_optional_number', '_positive', '_integer',
                            '_bounded_integer', '_tolerance', '_optional_number_or'})
INTEGER_READERS = frozenset({'_study_request_integer'})
OTHER_READERS = frozenset({'_boolean', '_flag', '_catalog_key', '_text', '_numbers',
                          '_items', '_objects', '_object'})


def dynamic_number_fields(path):
    feed = (set(transient._PRESSURANT_DEFAULTS) | set(transient._FEED_CONFIGURATION_DEFAULTS)
            | set(transient._FEED_STATE_DEFAULTS) | set(transient._ENGINE_CONDITION_DEFAULTS)
            | {'feed_line_resistance_pa_s2_kg2'})
    return {
        '/api/v1/transient/blowdown': (set(transient._BLOWDOWN_DEFAULTS), {'_blowdown_records'}),
        '/api/v1/transient/uncertainty': (set(transient._UNCERTAINTY_MODEL_DEFAULTS), {'_blowdown_records'}),
        '/api/v1/transient/regulated-feed': (feed, {'_feed_records'}),
        '/api/v1/transient/feed-line': (feed | set(transient._LINE_DEFAULTS),
                                       {'_feed_records', 'evaluate_feed_line_route'}),
        '/api/v1/transient/electric': (set(transient._ELECTRIC_NUMBERS), {'evaluate_electric_route'}),
        '/api/v1/transient/pulse-mode': (set(transient._PULSE_NUMBERS), {'plan_pulse_mode_route'}),
        '/api/v1/sass/mission': (set(transient._SASS_OVERRIDES), {'size_sass_mission_route'}),
    }.get(path, (set(), set()))


def numeric_payload_fields(path, root):
    """Return scalar and exact-integer fields, and any unresolved reader sites."""
    numbers, integers, seen, unresolved = set(), set(), set(), []
    dynamic, owners = dynamic_number_fields(path)

    def visit(function, parameter=None):
        function = inspect.unwrap(function)
        signature = inspect.signature(function)
        if parameter is None:
            parameter = next(iter(signature.parameters))
        marker = (id(function), parameter)
        if marker in seen:
            return
        seen.add(marker)
        tree = ast.parse(textwrap.dedent(inspect.getsource(function)))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Name):
                continue
            name = node.func.id
            keywords = {item.arg: item.value for item in node.keywords if item.arg is not None}
            positions = [i for i, arg in enumerate(node.args)
                         if isinstance(arg, ast.Name) and arg.id == parameter]
            named = [key for key, arg in keywords.items()
                     if isinstance(arg, ast.Name) and arg.id == parameter]
            if not positions and not named:
                continue
            if name in NUMBER_READERS | INTEGER_READERS and (0 in positions or 'payload' in named):
                target = integers if name in INTEGER_READERS else numbers
                key_node = node.args[1] if len(node.args) > 1 else keywords.get('key')
                if isinstance(key_node, ast.Constant) and isinstance(key_node.value, str):
                    target.add(key_node.value)
                elif (function.__name__ in owners and isinstance(key_node, ast.Name)
                      and key_node.id == 'key'):
                    target.update(dynamic)
                else:
                    unresolved.append((function.__name__, ast.unparse(node)))
            elif name not in OTHER_READERS:
                helper = function.__globals__.get(name)
                if inspect.isfunction(helper) and (
                        helper.__module__.startswith('rocket_propulsion.api.')
                        or helper.__module__ == root.__module__):
                    parameters = list(inspect.signature(helper).parameters)
                    for position in positions:
                        if position < len(parameters):
                            visit(helper, parameters[position])
                    for parameter_name in named:
                        if parameter_name in parameters:
                            visit(helper, parameter_name)

    visit(root)
    return numbers, integers, unresolved
