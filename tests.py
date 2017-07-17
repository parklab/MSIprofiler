from unittest import TestCase

from models import MicroSatelliteProfiler


class MSIProfilerTests(TestCase):
    def test_msi_profiler_no_args(self):
        with self.assertRaises(AttributeError):
            MicroSatelliteProfiler({})
