"""Prove that fault injection obeys controlled globals instead of false kills."""
from contextlib import ExitStack, contextmanager, redirect_stdout
import importlib
import importlib.util
import io
import json
import os
import sys
from threading import Thread
from pathlib import Path
import tempfile
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import patch

from scripts.mutation_runtime import compile_fault, evaluate_fault


SOURCE='''from functools import lru_cache
def hook(value):
    return True
def subject(value):
    if hook(value):
        return value+1
    return value
def recursive(depth):
    return recursive(depth-1)+1 if depth else 0
@lru_cache(maxsize=16)
def cached(value):
    if hook(value):
        return value+1
    return value
class Owner:
    @property
    def value(self):
        return 1 if hook(0) else 0
    @classmethod
    def read(cls,value):
        return value+1 if hook(value) else value
'''


def _temporary_modules(bindings):
    """Restore only owned fixture names; preserve imports completed elsewhere."""
    @contextmanager
    def manage():
        missing = object()
        before = {name: sys.modules.get(name, missing) for name in bindings}
        sys.modules.update(bindings)
        try:
            yield
        finally:
            for name, previous in before.items():
                if previous is missing:
                    sys.modules.pop(name, None)
                else:
                    sys.modules[name] = previous
    return manage()


class MutationRuntimeAttacks(unittest.TestCase):
    def test_module_fixture_cleanup_preserves_unrelated_concurrent_imports(self):
        owned = 'rocket_propulsion.studies.sass._owned_module_fixture'
        fresh = 'rocket_propulsion.studies.sass._new_module_fixture'
        late = 'rocket_propulsion.studies.sass._concurrently_loaded_fixture'
        sentinel = object()
        before = {name: sys.modules.get(name, sentinel) for name in (owned, fresh, late)}
        original, temporary, imported = (ModuleType(owned), ModuleType(owned), ModuleType(late))
        sys.modules[owned] = original
        sys.modules.pop(fresh, None)
        try:
            with _temporary_modules({owned: temporary, fresh: ModuleType(fresh)}):
                self.assertIs(sys.modules[owned], temporary)
                job = Thread(target=lambda: sys.modules.__setitem__(late, imported))
                job.start(); job.join(timeout=5)
                self.assertFalse(job.is_alive())
                self.assertIs(sys.modules[late], imported)
            self.assertIs(sys.modules[owned], original)
            self.assertNotIn(fresh, sys.modules)
            self.assertIs(sys.modules.get(late), imported,
                          'Fixture cleanup erased a dependency loaded by another thread.')
        finally:
            for name, previous in before.items():
                if previous is sentinel:
                    sys.modules.pop(name, None)
                else:
                    sys.modules[name] = previous

    def fixture(self,folder):
        name='rocket_propulsion.studies.sass._controlled_mutation_fixture'
        path=Path(folder)/'fixture.py';path.write_text(SOURCE)
        module=ModuleType(name);module.__file__=str(path)
        exec(compile(SOURCE,str(path),'exec'),module.__dict__)
        return module

    def exercise_runner(self,kind):
        if kind=='native' and importlib.util.find_spec('lunaris') is None:
            self.skipTest('Native mutation wiring requires the explicitly pinned native dependency context')
        runners={'feed':'scripts.run_feed_mutations','sass':'scripts.run_sass_mutations',
                 'native':'scripts.run_lunaris_impulse_mutations'}
        runner=importlib.import_module(runners[kind])
        original_import=importlib.import_module
        with tempfile.TemporaryDirectory() as folder:
            module=self.fixture(folder);original=module.subject;hook=module.hook
            (Path(folder)/'output/assurance').mkdir(parents=True)
            old_cwd=Path.cwd()
            try:
                os.chdir(folder)
                for active in (False,True):
                    class ControlledCase(unittest.TestCase):
                        def test_controlled(self):
                            with patch.object(module,'hook',return_value=active) as control:
                                self.assertEqual(module.subject(1),2 if active else 1)
                                control.assert_called_once_with(1)
                    baseline=unittest.TextTestRunner(stream=io.StringIO()).run(ControlledCase('test_controlled'))
                    self.assertTrue(baseline.wasSuccessful())
                    with ExitStack() as stack:
                        stack.enter_context(_temporary_modules({module.__name__:module}))
                        stack.enter_context(patch.object(runner,'ROOT',Path(folder)))
                        stack.enter_context(patch.object(runner,'verification_inputs',return_value={'fixture':'stable'}))
                        before,after='return value+1','return value+2'
                        if kind=='feed':
                            variants=(('controlled_fault',module,'subject',before,after,ControlledCase,'test_controlled'),)
                            output='feed_mutations.json'
                        elif kind=='sass':
                            variants=(('controlled_fault','fixture','subject',before,after,'test_controlled','adversarial'),)
                            def load(name):
                                if name=='rocket_propulsion.studies.sass.fixture':return module
                                if name=='tests.test_sass_adversarial':return SimpleNamespace(HostileRegressionTests=ControlledCase)
                                return original_import(name)
                            stack.enter_context(patch.object(runner.importlib,'import_module',side_effect=load))
                            output='mutations.json'
                        else:
                            variants=(('controlled_fault','subject',before,after,'test_controlled'),)
                            stack.enter_context(patch.object(runner,'plans',module))
                            stack.enter_context(patch.object(runner,'dependency_inputs',return_value={}))
                            stack.enter_context(patch.object(runner,'NativeLunarisImpulseAttacks',ControlledCase))
                            output='native_lunaris_mutations.json'
                        stack.enter_context(patch.object(runner,'MUTATIONS',variants))
                        with redirect_stdout(io.StringIO()):exit_status=runner.main()
                    row=json.loads((Path('output/assurance')/output).read_text())['mutations'][0]
                    self.assertEqual(row['errors'],0)
                    # With the controlled branch disabled, the injected fault
                    # cannot execute and must survive. A copied globals dict
                    # falsely executes the original hook and reports a kill.
                    self.assertEqual(row['status'],'killed' if active else 'survived')
                    self.assertEqual(exit_status,0 if active else 1)
                    self.assertIs(module.subject,original)
                    self.assertIs(module.hook,hook)
            finally:os.chdir(old_cwd)

    def test_feed_runner_respects_late_controlled_branch_and_restores_bindings(self):
        self.exercise_runner('feed')

    def test_sass_runner_respects_late_controlled_branch_and_restores_bindings(self):
        self.exercise_runner('sass')

    def test_native_runner_respects_late_controlled_branch_and_restores_bindings(self):
        self.exercise_runner('native')

    def test_recursive_fault_uses_the_same_live_binding_on_every_call(self):
        with tempfile.TemporaryDirectory() as folder:
            module=self.fixture(folder);original=module.recursive
            self.assertEqual(original(3),3)
            source,replacement=compile_fault(original,'+1 if','+2 if','<controlled-recursion>')
            self.assertIs(module.recursive,original)
            self.assertIs(replacement.__globals__,original.__globals__)
            with patch.object(module,'recursive',replacement):self.assertEqual(module.recursive(3),6)
            self.assertIs(module.recursive,original)

    def test_cached_function_and_property_do_not_capture_stale_callbacks(self):
        with tempfile.TemporaryDirectory() as folder:
            module=self.fixture(folder)
            _,cached=compile_fault(module.cached,'return value+1','return value+2','<controlled-cache>')
            _,property_=compile_fault(module.Owner.value.fget,'return 1 if','return 2 if','<controlled-property>')
            _,classmethod_=compile_fault(module.Owner.read,'value+1 if','value+2 if','<controlled-classmethod>')
            with patch.object(module,'hook',return_value=False):
                self.assertEqual(cached(1),1)
                self.assertEqual(property_.fget(module.Owner()),0)
                self.assertEqual(classmethod_.__func__(module.Owner,1),1)
            with patch.object(module,'hook',return_value=True):
                cached.cache_clear()
                self.assertEqual(cached(1),3)
                self.assertEqual(property_.fget(module.Owner()),2)
                self.assertEqual(classmethod_.__func__(module.Owner,1),3)

    def test_failing_original_control_cannot_be_counted_as_an_assertion_kill(self):
        with tempfile.TemporaryDirectory() as folder:
            module=self.fixture(folder);original=module.subject
            _,replacement=compile_fault(original,'return value+1','return value+2','<bad-control>')
            class WrongOracle(unittest.TestCase):
                def test_case(self):self.assertEqual(module.subject(1),100)
            bindings=unittest.mock.Mock(return_value=[(module,'subject')])
            report=evaluate_fault(bindings,original,replacement,WrongOracle,'test_case')
            self.assertEqual(report['status'],'invalid')
            self.assertFalse(report['baseline_passed'])
            self.assertFalse(report['execution_valid'])
            self.assertEqual(report['assertion_failures'],0)
            bindings.assert_not_called();self.assertIs(module.subject,original)

    def test_leaked_injected_state_fails_the_post_control_and_invalidates_the_kill(self):
        with tempfile.TemporaryDirectory() as folder:
            module=self.fixture(folder);original=module.subject;module.records=[]
            _,replacement=compile_fault(original,'return value+1',
                "records.append('dirty'); return value+2",'<state-leak>')
            class StatefulOracle(unittest.TestCase):
                def test_case(self):
                    self.assertEqual(module.records,[])
                    self.assertEqual(module.subject(1),2)
            report=evaluate_fault(lambda:[(module,'subject')],original,replacement,StatefulOracle,'test_case')
            self.assertTrue(report['baseline_passed'])
            self.assertFalse(report['post_control_passed'])
            self.assertEqual(report['assertion_failures'],1)
            self.assertEqual(report['status'],'invalid')
            self.assertIs(module.subject,original)

    def test_aliases_imported_by_controls_and_during_the_fault_are_patched_and_restored(self):
        with tempfile.TemporaryDirectory() as folder:
            module=self.fixture(folder);original=module.subject
            consumer=ModuleType('rocket_propulsion.studies.sass._late_control_consumer')
            new_consumer=ModuleType('rocket_propulsion.studies.sass._late_fault_consumer')
            module.new_consumer=new_consumer
            _,replacement=compile_fault(original,'return value+1',
                'new_consumer.alias=subject; return value+2','<late-alias>')
            class ImportingOracle(unittest.TestCase):
                def test_case(self):
                    if not hasattr(consumer,'alias'):consumer.alias=module.subject
                    self.assertEqual(consumer.alias(1),2)
            def bindings():
                return [(owner,key) for owner in (module,consumer) for key,value in list(vars(owner).items())
                        if value is original]
            with _temporary_modules({consumer.__name__:consumer,new_consumer.__name__:new_consumer}):
                report=evaluate_fault(bindings,original,replacement,ImportingOracle,'test_case')
            self.assertEqual(report['status'],'killed')
            self.assertTrue(report['baseline_passed'] and report['post_control_passed'])
            self.assertEqual(report['bindings_patched'],2)
            self.assertIs(consumer.alias,original);self.assertIs(new_consumer.alias,original)
