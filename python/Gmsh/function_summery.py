# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 13:05:32 2021

@author: benja
"""

# gives the node of ellements with dimention and the surfface name , from the geo file
def getNodes(dim=-1, tag=-1, includeBoundary=False, returnParametricCoord=True):
            """
            Get the nodes classified on the entity of dimension `dim' and tag `tag'. If
            `tag' < 0, get the nodes for all entities of dimension `dim'. If `dim' and
            `tag' are negative, get all the nodes in the mesh. `nodeTags' contains the
            node tags (their unique, strictly positive identification numbers). `coord'
            is a vector of length 3 times the length of `nodeTags' that contains the x,
            y, z coordinates of the nodes, concatenated: [n1x, n1y, n1z, n2x, ...]. If
            `dim' >= 0 and `returnParamtricCoord' is set, `parametricCoord' contains
            the parametric coordinates ([u1, u2, ...] or [u1, v1, u2, ...]) of the
            nodes, if available. The length of `parametricCoord' can be 0 or `dim'
            times the length of `nodeTags'. If `includeBoundary' is set, also return
            the nodes classified on the boundary of the entity (which will be
            reparametrized on the entity if `dim' >= 0 in order to compute their
            parametric coordinates).

            Return `nodeTags', `coord', `parametricCoord'.
            """
 # gives the node of ellements with type and the surfface name , from the geo file
def getNodesByElementType(elementType, tag=-1, returnParametricCoord=True):
            """

            Get the nodes classified on the entity of tag `tag', for all the elements
            of type `elementType'. The other arguments are treated as in `getNodes'.

            Return `nodeTags', `coord', `parametricCoord'.
            """

def getNode(nodeTag):
            """

            Get the coordinates and the parametric coordinates (if any) of the node
            with tag `tag'. This function relies on an internal cache (a vector in case
            of dense node numbering, a map otherwise); for large meshes accessing nodes
            in bulk is often preferable.

            Return `coord', `parametricCoord'.
            """

def getNodesForPhysicalGroup(dim, tag):
            """
            Get the nodes from all the elements belonging to the physical group of
            dimension `dim' and tag `tag'. `nodeTags' contains the node tags; `coord'
            is a vector of length 3 times the length of `nodeTags' that contains the x,
            y, z coordinates of the nodes, concatenated: [n1x, n1y, n1z, n2x, ...].

            Return `nodeTags', `coord'.
            """
def getElements(dim=-1, tag=-1):
            """
            Get the elements classified on the entity of dimension `dim' and tag `tag'.
            If `tag' < 0, get the elements for all entities of dimension `dim'. If
            `dim' and `tag' are negative, get all the elements in the mesh.
            `elementTypes' contains the MSH types of the elements (e.g. `2' for 3-node
            triangles: see `getElementProperties' to obtain the properties for a given
            element type). `elementTags' is a vector of the same length as
            `elementTypes'; each entry is a vector containing the tags (unique,
            strictly positive identifiers) of the elements of the corresponding type.
            `nodeTags' is also a vector of the same length as `elementTypes'; each
            entry is a vector of length equal to the number of elements of the given
            type times the number N of nodes for this type of element, that contains
            the node tags of all the elements of the given type, concatenated: [e1n1,
            e1n2, ..., e1nN, e2n1, ...].

            Return `elementTypes', `elementTags', `nodeTags'.
            """
def getElement(elementTag):
            """
            Get the type and node tags of the element with tag `tag'. This function
            relies on an internal cache (a vector in case of dense element numbering, a
            map otherwise); for large meshes accessing elements in bulk is often
            preferable.

            Return `elementType', `nodeTags'.
            """
 def getElementByCoordinates(x, y, z, dim=-1, strict=False):
            """
            Search the mesh for an element located at coordinates (`x', `y', `z'). This
            function performs a search in a spatial octree. If an element is found,
            return its tag, type and node tags, as well as the local coordinates (`u',
            `v', `w') within the reference element corresponding to search location. If
            `dim' is >= 0, only search for elements of the given dimension. If `strict'
            is not set, use a tolerance to find elements near the search location.

            Return `elementTag', `elementType', `nodeTags', `u', `v', `w'.
            """
def getLocalCoordinatesInElement(elementTag, x, y, z):
            """
            Return the local coordinates (`u', `v', `w') within the element
            `elementTag' corresponding to the model coordinates (`x', `y', `z'). This
            function relies on an internal cache (a vector in case of dense element
            numbering, a map otherwise); for large meshes accessing elements in bulk is
            often preferable.

            Return `u', `v', `w'.
            """
def getElementsByType(elementType, tag=-1, task=0, numTasks=1):
            """
            Get the elements of type `elementType' classified on the entity of tag
            `tag'. If `tag' < 0, get the elements for all entities. `elementTags' is a
            vector containing the tags (unique, strictly positive identifiers) of the
            elements of the corresponding type. `nodeTags' is a vector of length equal
            to the number of elements of the given type times the number N of nodes for
            this type of element, that contains the node tags of all the elements of
            the given type, concatenated: [e1n1, e1n2, ..., e1nN, e2n1, ...]. If
            `numTasks' > 1, only compute and return the part of the data indexed by
            `task'.

            Return `elementTags', `nodeTags'.
            """
def getEdges(nodeTags):
            """
            Get the global unique mesh edge identifiers `edgeTags' and orientations
            `edgeOrientation' for an input list of node tag pairs defining these edges,
            concatenated in the vector `nodeTags'.

            Return `edgeTags', `edgeOrientations'.
            """
def getFaces(faceType, nodeTags):
            """
            Get the global unique mesh face identifiers `faceTags' and orientations
            `faceOrientations' for an input list of node tag triplets (if `faceType' ==
            3) or quadruplets (if `faceType' == 4) defining these faces, concatenated
            in the vector `nodeTags'.

            Return `faceTags', `faceOrientations'.
            """
def getElementEdgeNodes(elementType, tag=-1, primary=False, task=0, numTasks=1):
            """
            Get the nodes on the edges of all elements of type `elementType' classified
            on the entity of tag `tag'. `nodeTags' contains the node tags of the edges
            for all the elements: [e1a1n1, e1a1n2, e1a2n1, ...]. Data is returned by
            element, with elements in the same order as in `getElements' and
            `getElementsByType'. If `primary' is set, only the primary (begin/end)
            nodes of the edges are returned. If `tag' < 0, get the edge nodes for all
            entities. If `numTasks' > 1, only compute and return the part of the data
            indexed by `task'.

            Return `nodeTags'.
            """
def getElementFaceNodes(elementType, faceType, tag=-1, primary=False, task=0, numTasks=1):
            """
            Get the nodes on the faces of type `faceType' (3 for triangular faces, 4
            for quadrangular faces) of all elements of type `elementType' classified on
            the entity of tag `tag'. `nodeTags' contains the node tags of the faces for
            all elements: [e1f1n1, ..., e1f1nFaceType, e1f2n1, ...]. Data is returned
            by element, with elements in the same order as in `getElements' and
            `getElementsByType'. If `primary' is set, only the primary (corner) nodes
            of the faces are returned. If `tag' < 0, get the face nodes for all
            entities. If `numTasks' > 1, only compute and return the part of the data
            indexed by `task'.

            Return `nodeTags'.
            """
