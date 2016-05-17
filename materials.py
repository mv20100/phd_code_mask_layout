class Materials:
	class Material:
		name = None
		resistivity = None

		def __init__(self,name,resistivity):
			self.name = name
			self.resistivity = resistivity

	gold = Material("Gold",2.44e-8)
	aluminium = Material("Aluminium",2.82e-8)
	copper = Material("Copper",1.68e-8)
	tungsten = Material("Tungsten",5.6e-8)
	platinum = Material("Platinum",1.06e-7)