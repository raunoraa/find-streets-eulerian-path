class Intersection:
    def __init__(self, id, movements, geometry, int_type):
        self.id = id
        self.movements = movements
        self.geometry = geometry
        self.int_type = int_type

    def __hash__(self):
        # Immutable identifier for hashing (id is unique)
        return hash(self.id)

    def __eq__(self, other):
        # Equality is based on ids
        return self.id == other.id

    def __repr__(self):
        return f"Intersection({self.id}, {self.movements}, {self.geometry}, {self.int_type})"
