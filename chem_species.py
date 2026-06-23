CHEM_ELEMENTS = ("H", "O", "Fe", "Mg", "C", "N", "Si", "S", "D")
CHEM_DESCRIPTOR_FIELDS = tuple(f"chem_{element}" for element in CHEM_ELEMENTS)
CHEM_STAR_FIELDS = tuple(f"{element}_star" for element in CHEM_ELEMENTS)
CHEM_GAS_FIELDS = tuple(f"{element}_gas" for element in CHEM_ELEMENTS)
