import _thread
import time

from iris_py import Iris_py

def print_time( threadName, delay):
   count = 0
   while count < 5:
      time.sleep(delay)
      count += 1
      print ("%s: %s" % ( threadName, time.ctime(time.time()) ))

try:
   _thread.start_new_thread(start_rx_stream, )
   _thread.start_new_thread(start_rx_stream, )
except:
   print ("Error: unable to start thread")
while 1:
   pass