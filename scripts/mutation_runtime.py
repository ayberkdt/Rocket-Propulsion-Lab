"""Compile one injected fault while preserving the callable's live module globals.

Definitions go to a separate locals dictionary. The function's global lookups
still resolve through its real module, so controlled fixtures, aliases and
recursive calls see the same later patches as the original implementation.
"""
import ast
import inspect
import textwrap
import sys
import io
import unittest
from contextlib import ExitStack
from unittest.mock import patch

MUTATION_RUNTIME_POLICY='live_module_globals_and_passing_before_after_controls_v1'


def compile_fault(function,before,after,label):
    source=textwrap.dedent(inspect.getsource(function))
    if source.count(before)!=1:
        raise RuntimeError(f'{label} must identify exactly one source expression.')
    original=inspect.unwrap(function)
    namespace=original.__globals__
    definitions={}
    exec(compile(ast.parse(source.replace(before,after)),label,'exec'),namespace,definitions)
    replacement=definitions[function.__name__]
    target=replacement.fget if isinstance(replacement,property) else getattr(replacement,'__func__',replacement)
    if inspect.unwrap(target).__globals__ is not namespace:
        raise RuntimeError(f'{label} detached the injected function from its live globals.')
    return source,replacement


def _modules():
    return [module for name,module in list(sys.modules.items()) if module is not None
            and name.startswith(('rocket_propulsion.','lunaris.','tests.','test_'))]


def _clear_caches():
    for module in _modules():
        for value in list(vars(module).values()):
            clear=getattr(value,'cache_clear',None)
            if callable(clear):clear()


def _run_case(test_class,test_name):
    _clear_caches()
    return unittest.TextTestRunner(stream=io.StringIO()).run(unittest.TestSuite((test_class(test_name),)))


def _single_execution(result):
    return result is not None and result.testsRun==1 and not (
        result.skipped or result.expectedFailures or result.unexpectedSuccesses)


def _control_passed(result):
    return _single_execution(result) and not result.failures and not result.errors


def evaluate_fault(bindings,original,replacement,test_class,test_name):
    """Require passing unmutated controls on both sides of one injected fault.

    Resolve aliases after the first control, which may import a new consumer.
    Aliases imported while the fault was active are restored before the second
    control as well. Control failures never become assertion kills.
    """
    before=_run_case(test_class,test_name)
    result=after=None
    bound=[]
    if _control_passed(before):
        bound=bindings()
        if not bound:raise RuntimeError('Injected fault has no live callable bindings.')
        with ExitStack() as stack:
            for owner,key in bound:stack.enter_context(patch.object(owner,key,replacement))
            result=_run_case(test_class,test_name)
        for module in _modules():
            for key,value in list(vars(module).items()):
                if value is replacement:setattr(module,key,original)
        after=_run_case(test_class,test_name)
    control_before,control_after=_control_passed(before),_control_passed(after)
    execution=_single_execution(result)
    failures=[] if result is None else result.failures
    errors=[] if result is None else result.errors
    status=('invalid' if not control_before or not control_after or not execution or errors
            else 'killed' if failures else 'survived')
    return {'status':status,'assertion_failures':len(failures),'errors':len(errors),
            'baseline_passed':control_before,'post_control_passed':control_after,
            'execution_valid':execution,'bindings_patched':len(bound),
            'failure_details':[detail for _,detail in failures],
            'error_details':[detail for _,detail in errors],
            'control_failure_details':[detail for r in (before,after) if r is not None
                                       for _,detail in r.failures+r.errors]}
