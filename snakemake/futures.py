
import os
import multiprocessing
import concurrent.futures
from concurrent.futures.process import _ResultItem


def _process_worker(call_queue, result_queue):
    """Evaluates calls from call_queue and places the results in result_queue.

    This worker is run in a separate process.

    Args:
        call_queue: A multiprocessing.Queue of _CallItems that will be read and
            evaluated by the worker.
        result_queue: A multiprocessing.Queue of _ResultItems that will written
            to by the worker.
        shutdown: A multiprocessing.Event that will be set as a signal to the
            worker that it should exit when call_queue is empty.
    """
    try:
        while True:
            call_item = call_queue.get(block=True)
            if call_item is None:
                # Wake up queue management thread
                result_queue.put(os.getpid())
                return
            try:
                r = call_item.fn(*call_item.args, **call_item.kwargs)
            except BaseException as e:
                result_queue.put(_ResultItem(call_item.work_id,
                                             exception=e))
            else:
                result_queue.put(_ResultItem(call_item.work_id,
                                             result=r))
    except KeyboardInterrupt:
        result_queue.put(None)


class ProcessPoolExecutor(concurrent.futures.ProcessPoolExecutor):
    def _adjust_process_count(self):
        for _ in range(len(self._processes), self._max_workers):
            p = multiprocessing.Process(
                    target=_process_worker,
                    args=(self._call_queue,
                          self._result_queue))
            p.start()
            self._processes[p.pid] = p
