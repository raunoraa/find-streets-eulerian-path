class Intersection:
    def __init__(self, id, movements, geometry):
        self.id = id
        self.movements = movements
        self.geometry = geometry
    
    def __hash__(self):
        # Immutable identifier for hashing (id is unique)
        return hash(self.id)
    
    def __eq__(self, other):
        # Equality is based on ids
        return self.id == other.id
    
    def __repr__(self):
        return f"Intersection({self.road_id}, {self.lane_id}, {self.is_entering})"