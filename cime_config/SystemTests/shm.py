"""
Implementation of the CAM SHM test.
This is a CAM specific test:
Verifies that MPI-3 node-shared memory for large read-only tables (cam_shmem_mod)
gives answers identical to per-rank private memory.
(1) Do a run with the default settings (node-shared tables if the shared-window
    probe passes; the atm log reports PASSED/FAILED)
(2) Do a run with cam_shmem_force_private = .true. (every rank holds private copies)
Where the probe fails (shared windows not usable), both runs use private memory.
"""
from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files
from CIME.status import append_testlog

logger = logging.getLogger(__name__)


class SHM(SystemTestsCompareTwo):
    def __init__(self, case, **kwargs):
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds = False,
                                       ignore_fieldlist_diffs = False,
                                       run_two_suffix = 'private',
                                       run_one_description = 'default cam_shmem (node-shared if available)',
                                       run_two_description = 'cam_shmem_force_private = .true.',
                                       **kwargs)

    def _case_one_setup(self):
        comments = "Run one: default cam_shmem settings (node-shared tables if the probe passes)."
        append_testlog(comments, self._orig_caseroot)

    def _case_two_setup(self):
        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "cam",
                                contents = "cam_shmem_force_private = .true.")
        comments = "Run two: cam_shmem_force_private = .true. (per-rank private memory)."
        append_testlog(comments, self._orig_caseroot)
