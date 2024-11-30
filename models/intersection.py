r"""
This submodule contains the Intersection class.
"""

class Intersection:
    """
    Represents an intersection with its unique ID, associated movements, geometry, and type.

    Attributes:
        id (int): The unique identifier of the intersection.
        movements (list): A list of movements (roads) associated with the intersection.
        geometry (shapely.geometry.base.BaseGeometry): A Shapely geometry object representing the intersection,
            which can be a Polygon, Point, or any geometry derived from Shapely's BaseGeometry class.
        int_type (str): The type of the intersection (will only be "intersection" in my implementation).

    Example:

    ::

        intersection = Intersection(
            id=0,
            movements=[
                    "Road #18005 -> Road #1",
                    "Road #18005 -> Road #11508"
                ],
            geometry=some_polygon,
            int_type="intersection"
        )
    """

    def __init__(self, id, movements, geometry, int_type):
        """
        Initialize an Intersection instance.

        Args:
            id (int): The unique identifier of the intersection.
            movements (list): A list of movements (roads).
            geometry (shapely.geometry.base.BaseGeometry): A Shapely geometry object representing the intersection,
                which can be a Polygon, Point, or any geometry derived from Shapely's BaseGeometry class.
            int_type (str): The type of the intersection (will only be "intersection" in my implementation).
        """
        self.id = id
        self.movements = movements
        self.geometry = geometry
        self.int_type = int_type

    def __hash__(self):
        """
        Generate a hash for the intersection based on its unique ID.

        Returns:
            int: The hash value of the intersection.

        Example:

        ::

            hash_value = hash(intersection)
        """
        return hash(self.id)

    def __eq__(self, other):
        """
        Check equality between two Intersection objects based on their IDs.

        Args:
            other (Intersection): Another Intersection instance to compare with.

        Returns:
            bool: True if both intersections have the same ID, False otherwise.

        Example:

        ::

            intersection1 == intersection2
        """
        return self.id == other.id

    def __repr__(self):
        """
        Return a string representation of the Intersection instance.

        Returns:
            str: A string describing the Intersection instance.

        Example:

        ::

            print(intersection)
            # Output: Intersection(0, ['Road #18005 -> Road #1', 'Road #18005 -> Road #11508'], <Polygon>, 'intersection')
        """
        return f"Intersection({self.id}, {self.movements}, {self.geometry}, {self.int_type})"
