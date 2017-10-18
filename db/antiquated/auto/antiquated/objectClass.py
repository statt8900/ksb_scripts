class Object:
	def toJSON(self): 
		import json
		return json.dumps(self, default=lambda o: o.__dict__,sort_keys=True, indent=4)
	def __eq__(self, other): #Equality of object means equality of all fields
		if type(other) is type(self):
			return self.__dict__ == other.__dict__
		return False
	def __ne__(self, other): return not self.__eq__(other)

