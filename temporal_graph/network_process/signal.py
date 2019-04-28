

__author__ = "Sumanta Mukherjee"
__all__ = ['SignalInput', 'StepSignal', 'PulseSignal', 'PeriodicPulse']


class SignalInput:
    def __init__(self, name='NoSignal'):
        self.__name = name

    def do_signal(self, t):
        return 0

    @property
    def name(self):
        return self.__name


class StepSignal(SignalInput):
    def __init__(self, t0, strength = 1):
        SignalInput.__init__(self, name='StepSignal')
        assert (t0 > 0)
        self._t0 = t0
        self._s = strength

    def do_signal(self, t):
        return self._s if self._t0 > t else 0


class PulseSignal(SignalInput):
    def __init__(self, ts, te, strength=1):
        SignalInput.__init__(self, name='PulseSignal')
        assert (ts < te) and (ts >= 0) and (te >=0)
        self._ts = ts
        self._te = te
        self._s = strength

    def do_signal(self, t):
        return self._s if (t > self._ts) and (t < self._te) else 0


class PeriodicPulse(SignalInput):
    def __init__(self, ts, span, periods, strength=1):
        SignalInput.__init__(self, name='PeriodicPulse')
        assert (span >= 0) and (periods > 0)
        self._ts = ts
        self._span = span
        self._period = periods
        self._s = strength

    def do_signal(self, t):
        dt = t - self._span
        dt = dt - int(dt / (self._span + self._period)) * (self._span + self._period)
        return self._s if dt < self._period else 0
