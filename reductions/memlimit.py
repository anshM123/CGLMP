"""memlimit.py -- put the current Windows process into a job object with a hard memory cap (default 1 GB).
Import and call cap() at the start of every script of this folder."""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import ctypes
from ctypes import wintypes


class IO_COUNTERS(ctypes.Structure):
    _fields_ = [(n, ctypes.c_ulonglong) for n in
                ("ReadOperationCount", "WriteOperationCount", "OtherOperationCount",
                 "ReadTransferCount", "WriteTransferCount", "OtherTransferCount")]


class JOBOBJECT_BASIC_LIMIT_INFORMATION(ctypes.Structure):
    _fields_ = [("PerProcessUserTimeLimit", ctypes.c_longlong),
                ("PerJobUserTimeLimit", ctypes.c_longlong),
                ("LimitFlags", wintypes.DWORD),
                ("MinimumWorkingSetSize", ctypes.c_size_t),
                ("MaximumWorkingSetSize", ctypes.c_size_t),
                ("ActiveProcessLimit", wintypes.DWORD),
                ("Affinity", ctypes.c_size_t),
                ("PriorityClass", wintypes.DWORD),
                ("SchedulingClass", wintypes.DWORD)]


class JOBOBJECT_EXTENDED_LIMIT_INFORMATION(ctypes.Structure):
    _fields_ = [("BasicLimitInformation", JOBOBJECT_BASIC_LIMIT_INFORMATION),
                ("IoInfo", IO_COUNTERS),
                ("ProcessMemoryLimit", ctypes.c_size_t),
                ("JobMemoryLimit", ctypes.c_size_t),
                ("PeakProcessMemoryUsed", ctypes.c_size_t),
                ("PeakJobMemoryUsed", ctypes.c_size_t)]


_job = None


def cap(mb=1000):
    global _job
    try:
        k32 = ctypes.windll.kernel32
        k32.CreateJobObjectW.restype = wintypes.HANDLE
        job = k32.CreateJobObjectW(None, None)
        info = JOBOBJECT_EXTENDED_LIMIT_INFORMATION()
        JOB_OBJECT_LIMIT_PROCESS_MEMORY = 0x00000100
        info.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_PROCESS_MEMORY
        info.ProcessMemoryLimit = mb * 1024 * 1024
        ok = k32.SetInformationJobObject(wintypes.HANDLE(job), 9, ctypes.byref(info), ctypes.sizeof(info))
        hproc = k32.GetCurrentProcess()
        ok2 = k32.AssignProcessToJobObject(wintypes.HANDLE(job), wintypes.HANDLE(hproc))
        _job = job
        return bool(ok and ok2)
    except Exception:
        return False


if __name__ == "__main__":
    print("cap:", cap(200))
    import numpy as np
    try:
        a = np.ones((400, 1024, 1024 // 8))  # ~400 MB
        print("allocated (unexpected)")
    except MemoryError:
        print("MemoryError as expected")
