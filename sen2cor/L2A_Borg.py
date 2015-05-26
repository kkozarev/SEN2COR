'''
Created on Feb 22, 2012
The Borg Pattern for sharing states
@author: umuellerwilm
'''
class Borg(object):
   _shared = {}
   def __new__(cls, *p, **k):
      inst = object.__new__(cls)
      inst.__dict__ = cls._shared
      return inst

