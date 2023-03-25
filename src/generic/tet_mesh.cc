// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================

#include <algorithm>

#include "map_matrix.h"
#include "tet_mesh.h"
#include "Telements.h"


namespace oomph
{
  //=======================================================================
  /// Constructor for a FacetedSurface created from a list of nodes
  /// and connectivity information. This is used in remeshing
  //=======================================================================
  TetMeshFacetedClosedSurfaceForRemesh::TetMeshFacetedClosedSurfaceForRemesh(
    Vector<Node*> const& vertex_node_pt,
    Vector<Vector<unsigned>> const& facet_connectivity,
    Vector<unsigned> const& facet_boundary_id)
    : TetMeshFacetedClosedSurface()
  {
    // Create the vertices
    unsigned n_vertex = vertex_node_pt.size();
    Vertex_pt.resize(n_vertex);
    for (unsigned v = 0; v < n_vertex; ++v)
    {
      Vertex_pt[v] = new TetMeshVertex(vertex_node_pt[v]);
    }

    // Create the facets
    unsigned n_facet = facet_connectivity.size();
    Facet_pt.resize(n_facet);
    for (unsigned f = 0; f < n_facet; ++f)
    {
      unsigned n_vertex_on_facet = facet_connectivity[f].size();
      Facet_pt[f] = new TetMeshFacet(n_vertex_on_facet);
      for (unsigned i = 0; i < n_vertex_on_facet; ++i)
      {
        Facet_pt[f]->set_vertex_pt(i, Vertex_pt[facet_connectivity[f][i]]);
      }
      // Add in the boundary id
      Facet_pt[f]->set_one_based_boundary_id(facet_boundary_id[f]);
    }
  }

  //=================================================================
  /// Destructor. Delete allocated memory
  //================================================================
  TetMeshFacetedClosedSurfaceForRemesh::~TetMeshFacetedClosedSurfaceForRemesh()
  {
    // Delete the facets and the vertices
    unsigned n_facet = this->nfacet();
    for (unsigned f = 0; f < n_facet; f++)
    {
      delete Facet_pt[f];
    }
    unsigned n_vertex = this->nvertex();
    for (unsigned v = 0; v < n_vertex; v++)
    {
      delete Vertex_pt[v];
    }
  }


  //================================================================
  /// Global static data that specifies the permitted
  /// error in the setup of the boundary coordinates
  //================================================================
  double TetMeshBase::Tolerance_for_boundary_finding = 1.0e-5;


  //================================================================
  /// Setup lookup schemes which establish which elements are located
  /// next to which boundaries (Doc to outfile if it's open).
  //================================================================
  void TetMeshBase::setup_boundary_element_info(std::ostream& outfile)
  {
    // Should we document the output here
    bool doc = false;

    if (outfile) doc = true;

    // Number of boundaries
    unsigned nbound = nboundary();

    // Wipe/allocate storage for arrays
    Boundary_element_pt.clear();
    Face_index_at_boundary.clear();
    Boundary_element_pt.resize(nbound);
    Face_index_at_boundary.resize(nbound);
    Lookup_for_elements_next_boundary_is_setup = false;

    // Temporary vector of vectors of pointers to elements on the boundaries:
    // This is a vector to ensure UNIQUE ordering in all processors
    Vector<Vector<FiniteElement*>> vector_of_boundary_element_pt;
    vector_of_boundary_element_pt.resize(nbound);

    // Matrix map for working out the fixed face for elements on boundary
    MapMatrixMixed<unsigned, FiniteElement*, int> face_identifier;

    // Loop over elements
    //-------------------
    unsigned nel = nelement();


    // Get pointer to vector of boundaries that the
    // node lives on
    Vector<std::set<unsigned>*> boundaries_pt(4, 0);

    for (unsigned e = 0; e < nel; e++)
    {
      // Get pointer to element
      FiniteElement* fe_pt = finite_element_pt(e);

      if (doc) outfile << "Element: " << e << " " << fe_pt << std::endl;

      // Only include 3D elements! Some meshes contain interface elements too.
      if (fe_pt->dim() == 3)
      {
        // Loop over the element's nodes and find out which boundaries they're
        // on
        // ----------------------------------------------------------------------
        // We need only loop over the corner nodes
        for (unsigned i = 0; i < 4; i++)
        {
          fe_pt->node_pt(i)->get_boundaries_pt(boundaries_pt[i]);
        }

        // Find the common boundaries of each face
        Vector<std::set<unsigned>> face(4);

        // NOTE: Face indices defined in Telements.h

        // Face 3 connnects points 0, 1 and 2
        if (boundaries_pt[0] && boundaries_pt[1] && boundaries_pt[2])
        {
          std::set<unsigned> aux;

          std::set_intersection(
            boundaries_pt[0]->begin(),
            boundaries_pt[0]->end(),
            boundaries_pt[1]->begin(),
            boundaries_pt[1]->end(),
            std::insert_iterator<std::set<unsigned>>(aux, aux.begin()));

          std::set_intersection(
            aux.begin(),
            aux.end(),
            boundaries_pt[2]->begin(),
            boundaries_pt[2]->end(),
            std::insert_iterator<std::set<unsigned>>(face[3], face[3].begin()));
        }

        // Face 2 connects points 0, 1 and 3
        if (boundaries_pt[0] && boundaries_pt[1] && boundaries_pt[3])
        {
          std::set<unsigned> aux;

          std::set_intersection(
            boundaries_pt[0]->begin(),
            boundaries_pt[0]->end(),
            boundaries_pt[1]->begin(),
            boundaries_pt[1]->end(),
            std::insert_iterator<std::set<unsigned>>(aux, aux.begin()));

          std::set_intersection(
            aux.begin(),
            aux.end(),
            boundaries_pt[3]->begin(),
            boundaries_pt[3]->end(),
            std::insert_iterator<std::set<unsigned>>(face[2], face[2].begin()));
        }

        // Face 1 connects points 0, 2 and 3
        if (boundaries_pt[0] && boundaries_pt[2] && boundaries_pt[3])
        {
          std::set<unsigned> aux;

          std::set_intersection(
            boundaries_pt[0]->begin(),
            boundaries_pt[0]->end(),
            boundaries_pt[2]->begin(),
            boundaries_pt[2]->end(),
            std::insert_iterator<std::set<unsigned>>(aux, aux.begin()));

          std::set_intersection(
            aux.begin(),
            aux.end(),
            boundaries_pt[3]->begin(),
            boundaries_pt[3]->end(),
            std::insert_iterator<std::set<unsigned>>(face[1], face[1].begin()));
        }

        // Face 0 connects points 1, 2 and 3
        if (boundaries_pt[1] && boundaries_pt[2] && boundaries_pt[3])
        {
          std::set<unsigned> aux;

          std::set_intersection(
            boundaries_pt[1]->begin(),
            boundaries_pt[1]->end(),
            boundaries_pt[2]->begin(),
            boundaries_pt[2]->end(),
            std::insert_iterator<std::set<unsigned>>(aux, aux.begin()));

          std::set_intersection(
            aux.begin(),
            aux.end(),
            boundaries_pt[3]->begin(),
            boundaries_pt[3]->end(),
            std::insert_iterator<std::set<unsigned>>(face[0], face[0].begin()));
        }


        // We now know whether any faces lay on the boundaries
        for (unsigned i = 0; i < 4; i++)
        {
          // How many boundaries are there
          unsigned count = 0;

          // The number of the boundary
          int boundary = -1;

          // Loop over all the members of the set and add to the count
          // and set the boundary
          for (std::set<unsigned>::iterator it = face[i].begin();
               it != face[i].end();
               ++it)
          {
            ++count;
            boundary = *it;
          }

          // If we're on more than one boundary, this is weird, so die
          if (count > 1)
          {
            std::ostringstream error_stream;
            fe_pt->output(error_stream);
            error_stream << "Face " << i << " is on " << count
                         << " boundaries.\n";
            error_stream << "This is rather strange.\n";
            error_stream << "Your mesh may be too coarse or your tet mesh\n";
            error_stream << "may be screwed up. I'm skipping the automated\n";
            error_stream << "setup of the elements next to the boundaries\n";
            error_stream << "lookup schemes.\n";
            OomphLibWarning(error_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
          }

          // If we have a boundary then add this to the appropriate set
          if (boundary >= 0)
          {
            // Does the pointer already exits in the vector
            Vector<FiniteElement*>::iterator b_el_it = std::find(
              vector_of_boundary_element_pt[static_cast<unsigned>(boundary)]
                .begin(),
              vector_of_boundary_element_pt[static_cast<unsigned>(boundary)]
                .end(),
              fe_pt);

            // Only insert if we have not found it (i.e. got to the end)
            if (b_el_it ==
                vector_of_boundary_element_pt[static_cast<unsigned>(boundary)]
                  .end())
            {
              vector_of_boundary_element_pt[static_cast<unsigned>(boundary)]
                .push_back(fe_pt);
            }

            // Also set the fixed face
            face_identifier(static_cast<unsigned>(boundary), fe_pt) = i;
          }
        }

        // Now we set the pointers to the boundary sets to zero
        for (unsigned i = 0; i < 4; i++)
        {
          boundaries_pt[i] = 0;
        }
      }
    }

    // Now copy everything across into permanent arrays
    //-------------------------------------------------

    // Loop over boundaries
    //---------------------
    for (unsigned i = 0; i < nbound; i++)
    {
      // Number of elements on this boundary (currently stored in a set)
      unsigned nel = vector_of_boundary_element_pt[i].size();

      // Allocate storage for the coordinate identifiers
      Face_index_at_boundary[i].resize(nel);

      unsigned e_count = 0;
      typedef Vector<FiniteElement*>::iterator IT;
      for (IT it = vector_of_boundary_element_pt[i].begin();
           it != vector_of_boundary_element_pt[i].end();
           it++)
      {
        // Recover pointer to element
        FiniteElement* fe_pt = *it;

        // Add to permanent storage
        Boundary_element_pt[i].push_back(fe_pt);

        Face_index_at_boundary[i][e_count] = face_identifier(i, fe_pt);

        // Increment counter
        e_count++;
      }
    }


    // Doc?
    //-----
    if (doc)
    {
      // Loop over boundaries
      for (unsigned i = 0; i < nbound; i++)
      {
        unsigned nel = Boundary_element_pt[i].size();
        outfile << "Boundary: " << i << " is adjacent to " << nel << " elements"
                << std::endl;

        // Loop over elements on given boundary
        for (unsigned e = 0; e < nel; e++)
        {
          FiniteElement* fe_pt = Boundary_element_pt[i][e];
          outfile << "Boundary element:" << fe_pt
                  << " Face index of boundary is "
                  << Face_index_at_boundary[i][e] << std::endl;
        }
      }
    }


    // Lookup scheme has now been setup yet
    Lookup_for_elements_next_boundary_is_setup = true;
  }


  //======================================================================
  /// Assess mesh quality: Ratio of max. edge length to min. height,
  /// so if it's very large it's BAAAAAD.
  //======================================================================
  void TetMeshBase::assess_mesh_quality(std::ofstream& some_file)
  {
    Vector<Vector<double>> edge(6);
    for (unsigned e = 0; e < 6; e++)
    {
      edge[e].resize(3);
    }
    unsigned nel = this->nelement();
    for (unsigned e = 0; e < nel; e++)
    {
      FiniteElement* fe_pt = this->finite_element_pt(e);
      for (unsigned i = 0; i < 3; i++)
      {
        edge[0][i] = fe_pt->node_pt(2)->x(i) - fe_pt->node_pt(1)->x(i);
        edge[1][i] = fe_pt->node_pt(0)->x(i) - fe_pt->node_pt(2)->x(i);
        edge[2][i] = fe_pt->node_pt(1)->x(i) - fe_pt->node_pt(0)->x(i);
        edge[3][i] = fe_pt->node_pt(3)->x(i) - fe_pt->node_pt(0)->x(i);
        edge[4][i] = fe_pt->node_pt(3)->x(i) - fe_pt->node_pt(1)->x(i);
        edge[5][i] = fe_pt->node_pt(3)->x(i) - fe_pt->node_pt(2)->x(i);
      }

      double max_length = 0.0;
      for (unsigned j = 0; j < 6; j++)
      {
        double length = 0.0;
        for (unsigned i = 0; i < 3; i++)
        {
          length += edge[j][i] * edge[j][i];
        }
        length = sqrt(length);
        if (length > max_length) max_length = length;
      }


      double min_height = DBL_MAX;
      for (unsigned j = 0; j < 4; j++)
      {
        Vector<double> normal(3);
        unsigned e0 = 0;
        unsigned e1 = 0;
        unsigned e2 = 0;
        switch (j)
        {
          case 0:
            e0 = 4;
            e1 = 5;
            e2 = 1;
            break;

          case 1:
            e0 = 1;
            e1 = 3;
            e2 = 2;
            break;

          case 2:
            e0 = 3;
            e1 = 4;
            e2 = 1;
            break;

          case 3:
            e0 = 1;
            e1 = 2;
            e2 = 3;
            break;

          default:

            oomph_info << "never get here\n";
            abort();
        }

        normal[0] = edge[e0][1] * edge[e1][2] - edge[e0][2] * edge[e1][1];
        normal[1] = edge[e0][2] * edge[e1][0] - edge[e0][0] * edge[e1][2];
        normal[2] = edge[e0][0] * edge[e1][1] - edge[e0][1] * edge[e1][0];
        double norm =
          normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
        double inv_norm = 1.0 / sqrt(norm);
        normal[0] *= inv_norm;
        normal[1] *= inv_norm;
        normal[2] *= inv_norm;

        double height = fabs(edge[e2][0] * normal[0] + edge[e2][1] * normal[1] +
                             edge[e2][2] * normal[2]);

        if (height < min_height) min_height = height;
      }

      double aspect_ratio = max_length / min_height;

      some_file << "ZONE N=4, E=1, F=FEPOINT, ET=TETRAHEDRON\n";
      for (unsigned j = 0; j < 4; j++)
      {
        for (unsigned i = 0; i < 3; i++)
        {
          some_file << fe_pt->node_pt(j)->x(i) << " ";
        }
        some_file << aspect_ratio << std::endl;
      }
      some_file << "1 2 3 4" << std::endl;
    }
  }


  //======================================================================
  /// Move the nodes on boundaries with associated Geometric Objects (if any)
  /// so that they exactly coincide with the geometric object. This requires
  /// that the boundary coordinates are set up consistently
  //======================================================================
  void TetMeshBase::snap_nodes_onto_geometric_objects()
  {
    // Backup in case elements get inverted
    std::map<Node*, Vector<double>> old_nodal_posn;
    std::map<Node*, Vector<double>> new_nodal_posn;
    unsigned nnod = nnode();
    for (unsigned j = 0; j < nnod; j++)
    {
      Node* nod_pt = node_pt(j);
      Vector<double> x(3);
      nod_pt->position(x);
      old_nodal_posn[nod_pt] = x;
    }

    // Loop over all boundaries
    unsigned n_bound = this->nboundary();
    for (unsigned b = 0; b < n_bound; b++)
    {
      bool do_it = true;

      // Accumulate reason for not snapping
      std::stringstream reason;
      reason << "Can't snap nodes on boundary " << b
             << " onto geom object because: \n";

      TetMeshFacetedSurface* faceted_surface_pt = 0;
      std::map<unsigned, TetMeshFacetedSurface*>::iterator it =
        Tet_mesh_faceted_surface_pt.find(b);
      if (it != Tet_mesh_faceted_surface_pt.end())
      {
        faceted_surface_pt = (*it).second;
      }

      // Facet associated with this boundary?
      if (faceted_surface_pt == 0)
      {
        reason << "-- no facets asssociated with boundary\n";
        do_it = false;
      }

      // Get geom object associated with facet
      GeomObject* geom_obj_pt = 0;
      if (do_it)
      {
        geom_obj_pt = faceted_surface_pt->geom_object_with_boundaries_pt();
        if (geom_obj_pt == 0)
        {
          reason << "-- no geom object associated with boundary\n";
          do_it = false;
        }
      }

      // Triangular facet?
      if (Triangular_facet_vertex_boundary_coordinate[b].size() == 0)
      {
        reason << "-- facet has to be triangular and vertex coordinates have\n"
               << "   to have been set up\n";
        do_it = false;
      }

      // We need boundary coordinates!
      if (!Boundary_coordinate_exists[b])
      {
        reason << "-- no boundary coordinates were set up\n";
        do_it = false;
      }


      // Which facet is associated with this boundary?
      unsigned facet_id_of_boundary = 0;
      TetMeshFacet* f_pt = 0;
      if (do_it)
      {
        unsigned nf = faceted_surface_pt->nfacet();
        for (unsigned f = 0; f < nf; f++)
        {
          if ((faceted_surface_pt->one_based_facet_boundary_id(f) - 1) == b)
          {
            facet_id_of_boundary = f;
            break;
          }
        }
        f_pt = faceted_surface_pt->facet_pt(facet_id_of_boundary);


        // Three vertices?
        unsigned nv = f_pt->nvertex();
        if (nv != 3)
        {
          reason << "-- number of facet vertices is " << nv
                 << " rather than 3\n";
          do_it = false;
        }

        // Have we set up zeta coordinates in geometric object?
        if ((f_pt->vertex_pt(0)->zeta_in_geom_object().size() != 2) ||
            (f_pt->vertex_pt(1)->zeta_in_geom_object().size() != 2) ||
            (f_pt->vertex_pt(2)->zeta_in_geom_object().size() != 2))
        {
          reason << "-- no boundary coordinates were set up\n";
          do_it = false;
        }
      }


      // Are we ready to go?
      if (!do_it)
      {
        const bool tell_us_why = false;
        if (tell_us_why)
        {
          oomph_info << reason.str() << std::endl;
        }
      }
      else
      {
        // Setup area coordinantes in triangular facet
        double x1 = Triangular_facet_vertex_boundary_coordinate[b][0][0];
        double y1 = Triangular_facet_vertex_boundary_coordinate[b][0][1];

        double x2 = Triangular_facet_vertex_boundary_coordinate[b][1][0];
        double y2 = Triangular_facet_vertex_boundary_coordinate[b][1][1];

        double x3 = Triangular_facet_vertex_boundary_coordinate[b][2][0];
        double y3 = Triangular_facet_vertex_boundary_coordinate[b][2][1];

        double detT = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);


        // Boundary coordinate (cartesian coordinates inside facet)
        Vector<double> zeta(2);

        // Loop over all nodes on that boundary
        const unsigned n_boundary_node = this->nboundary_node(b);
        for (unsigned n = 0; n < n_boundary_node; ++n)
        {
          // Get the boundary node and coordinates
          Node* const nod_pt = this->boundary_node_pt(b, n);
          nod_pt->get_coordinates_on_boundary(b, zeta);

          // Now we have zeta, the cartesian boundary coordinates
          // in the (assumed to be triangular) boundary facet; let's
          // work out the area coordinates
          // Notation as in
          // https://en.wikipedia.org/wiki/Barycentric_coordinate_system
          double s0 =
            ((y2 - y3) * (zeta[0] - x3) + (x3 - x2) * (zeta[1] - y3)) / detT;
          double s1 =
            ((y3 - y1) * (zeta[0] - x3) + (x1 - x3) * (zeta[1] - y3)) / detT;
          double s2 = 1.0 - s0 - s1;

          Vector<double> zeta_in_geom_obj(2, 0.0);
          Vector<double> position_from_geom_obj(3, 0.0);

          // Vertex zeta coordinates
          Vector<double> zeta_0(2);
          zeta_0 = f_pt->vertex_pt(0)->zeta_in_geom_object();

          Vector<double> zeta_1(2);
          zeta_1 = f_pt->vertex_pt(1)->zeta_in_geom_object();

          Vector<double> zeta_2(2);
          zeta_2 = f_pt->vertex_pt(2)->zeta_in_geom_object();


#ifdef PARANOID

          // Compute zeta values of the vertices from parametrisation of
          // boundaries
          double tol = 1.0e-12;
          Vector<double> zeta_from_boundary(2);
          Vector<double> zeta_vertex(2);
          for (unsigned v = 0; v < 3; v++)
          {
            zeta_vertex = f_pt->vertex_pt(v)->zeta_in_geom_object();
            for (unsigned alt = 0; alt < 2; alt++)
            {
              switch (v)
              {
                case 0:

                  if (alt == 0)
                  {
                    faceted_surface_pt->boundary_zeta01(
                      facet_id_of_boundary, 0.0, zeta_from_boundary);
                  }
                  else
                  {
                    faceted_surface_pt->boundary_zeta20(
                      facet_id_of_boundary, 1.0, zeta_from_boundary);
                  }
                  break;

                case 1:

                  if (alt == 0)
                  {
                    faceted_surface_pt->boundary_zeta01(
                      facet_id_of_boundary, 1.0, zeta_from_boundary);
                  }
                  else
                  {
                    faceted_surface_pt->boundary_zeta12(
                      facet_id_of_boundary, 0.0, zeta_from_boundary);
                  }
                  break;

                case 2:

                  if (alt == 0)
                  {
                    faceted_surface_pt->boundary_zeta12(
                      facet_id_of_boundary, 1.0, zeta_from_boundary);
                  }
                  else
                  {
                    faceted_surface_pt->boundary_zeta20(
                      facet_id_of_boundary, 0.0, zeta_from_boundary);
                  }
                  break;
              }

              double error =
                sqrt(pow((zeta_vertex[0] - zeta_from_boundary[0]), 2) +
                     pow((zeta_vertex[1] - zeta_from_boundary[1]), 2));
              if (error > tol)
              {
                std::ostringstream error_message;
                error_message
                  << "Error in parametrisation of boundary coordinates \n"
                  << "for vertex " << v << " [alt=" << alt << "] in facet "
                  << facet_id_of_boundary << " : \n"
                  << "zeta_vertex = [ " << zeta_vertex[0] << " "
                  << zeta_vertex[1] << " ] \n"
                  << "zeta_from_boundary      = [ " << zeta_from_boundary[0]
                  << " " << zeta_from_boundary[1] << " ] \n"
                  << std::endl;
                throw OomphLibError(error_message.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
            }
          }

#endif

          // Compute zeta values of the interpolation parameters
          Vector<double> zeta_a(3, 0.0);
          Vector<double> zeta_b(3, 0.0);
          Vector<double> zeta_c(3, 0.0);
          Vector<double> zeta_d(3, 0.0);
          Vector<double> zeta_e(3, 0.0);
          Vector<double> zeta_f(3, 0.0);
          faceted_surface_pt->boundary_zeta01(facet_id_of_boundary, s1, zeta_a);
          faceted_surface_pt->boundary_zeta01(
            facet_id_of_boundary, 1.0 - s0, zeta_d);

          faceted_surface_pt->boundary_zeta12(facet_id_of_boundary, s2, zeta_c);
          faceted_surface_pt->boundary_zeta12(
            facet_id_of_boundary, 1.0 - s1, zeta_f);

          faceted_surface_pt->boundary_zeta20(
            facet_id_of_boundary, 1.0 - s2, zeta_b);
          faceted_surface_pt->boundary_zeta20(facet_id_of_boundary, s0, zeta_e);

          // Transfinite mapping
          zeta_in_geom_obj[0] = s0 * (zeta_a[0] + zeta_b[0] - zeta_0[0]) +
                                s1 * (zeta_c[0] + zeta_d[0] - zeta_1[0]) +
                                s2 * (zeta_e[0] + zeta_f[0] - zeta_2[0]);
          zeta_in_geom_obj[1] = s0 * (zeta_a[1] + zeta_b[1] - zeta_0[1]) +
                                s1 * (zeta_c[1] + zeta_d[1] - zeta_1[1]) +
                                s2 * (zeta_e[1] + zeta_f[1] - zeta_2[1]);

          unsigned n_tvalues =
            1 + nod_pt->position_time_stepper_pt()->nprev_values();
          for (unsigned t = 0; t < n_tvalues; ++t)
          {
            // Get the position according to the underlying geometric object
            geom_obj_pt->position(t, zeta_in_geom_obj, position_from_geom_obj);

            // Move the node
            for (unsigned i = 0; i < 3; i++)
            {
              nod_pt->x(t, i) = position_from_geom_obj[i];
            }
          }
        }
      }
    }

    // Check if any element is inverted
    bool some_element_is_inverted = false;
    unsigned count = 0;
    unsigned nel = nelement();
    for (unsigned e = 0; e < nel; e++)
    {
      FiniteElement* el_pt = finite_element_pt(e);
      bool passed = true;
      el_pt->check_J_eulerian_at_knots(passed);
      if (!passed)
      {
        some_element_is_inverted = true;
        char filename[100];
        std::ofstream some_file;
        sprintf(filename, "overly_distorted_element%i.dat", count);
        some_file.open(filename);
        unsigned nnod_1d = el_pt->nnode_1d();
        el_pt->output(some_file, nnod_1d);
        some_file.close();

        // Reset to old nodal position
        unsigned nnod = el_pt->nnode();
        for (unsigned j = 0; j < nnod; j++)
        {
          Node* nod_pt = el_pt->node_pt(j);
          Vector<double> x_current(3);
          nod_pt->position(x_current);
          new_nodal_posn[nod_pt] = x_current;
          Vector<double> old_x(old_nodal_posn[nod_pt]);
          for (unsigned i = 0; i < 3; i++)
          {
            nod_pt->x(i) = old_x[i];
          }
        }

        // Plot
        sprintf(filename, "orig_overly_distorted_element%i.dat", count);
        some_file.open(filename);
        el_pt->output(some_file, nnod_1d);
        some_file.close();

        // Reset
        for (unsigned j = 0; j < nnod; j++)
        {
          Node* nod_pt = el_pt->node_pt(j);
          for (unsigned i = 0; i < 3; i++)
          {
            nod_pt->x(i) = new_nodal_posn[nod_pt][i];
          }
        }

        // Bump
        count++;
      }
    }
    if (some_element_is_inverted)
    {
      std::ostringstream error_message;
      error_message
        << "A number of elements, namely: " << count
        << " are inverted after snapping. Their shapes are in "
        << " overly_distorted_element*.dat and "
           "orig_overly_distorted_element*.dat"
        << "Next person to get this error: Please implement a straightforward\n"
        << "variant of one of the functors in src/mesh_smoothing to switch\n"
        << "to harmonic mapping\n"
        << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    else
    {
      oomph_info << "No elements are inverted after snapping. Yay!"
                 << std::endl;
    }
  }


  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////

  // A series of examples of TetMeshFacetedSurfaces

  //=============================================================================
  /// CubicTetMeshFacetedClosedSurface from (+/- half width)^2 x +/- half length
  //=============================================================================
  class CubicTetMeshFacetedSurface : public virtual TetMeshFacetedClosedSurface
  {
  public:
    /// Constructor: The six faces are enumerated consecutively
    /// starting from boundary_id_offset (defaults to zero in which
    /// case the boundaries are enumerated (in one-based fashion) 1,...,6.
    CubicTetMeshFacetedSurface(const double& box_half_width,
                               const double& box_half_length,
                               const unsigned& one_based_boundary_id_offset = 0)
    {
      unsigned one_based_region_id = 0;
      Vector<double> offset(3, 0.0);
      build_it(one_based_region_id,
               box_half_width,
               box_half_length,
               offset,
               one_based_boundary_id_offset);
    }

    /// Constructor: The six faces are enumerated consecutively
    /// starting from boundary_id_offset (defaults to zero in which
    /// case the boundaries are enumerated (in one-based fashion) 1,...,6.
    CubicTetMeshFacetedSurface(const double& box_half_width,
                               const double& box_half_length,
                               const Vector<double>& offset,
                               const unsigned& one_based_boundary_id_offset = 0)
    {
      unsigned one_based_region_id = 0;
      build_it(one_based_region_id,
               box_half_width,
               box_half_length,
               offset,
               one_based_boundary_id_offset);
    }
    /// Constructor: Faceted surface defines region with (one-based!)
    /// ID one_based_region_id. If this parameter is zero, no region ID is
    /// assigned. The six faces are enumerated consecutively
    /// starting from boundary_id_offset (defaults to zero in which
    /// case the boundaries are enumerated (in one-based fashion) 1,...,6.
    CubicTetMeshFacetedSurface(const unsigned& one_based_region_id,
                               const double& box_half_width,
                               const double& box_half_length,
                               const unsigned& one_based_boundary_id_offset = 0)
    {
      Vector<double> offset(3, 0.0);
      build_it(one_based_region_id,
               box_half_width,
               box_half_length,
               offset,
               one_based_boundary_id_offset);
    }

    /// Constructor: Faceted surface defines region with (one-based!)
    /// ID one_based_region_id. If this parameter is zero, no region ID is
    /// assigned. The six faces are enumerated consecutively
    /// starting from boundary_id_offset (defaults to zero in which
    /// case the boundaries are enumerated (in one-based fashion) 1,...,6.
    CubicTetMeshFacetedSurface(const unsigned& one_based_region_id,
                               const double& box_half_width,
                               const double& box_half_length,
                               const Vector<double>& offset,
                               const unsigned& one_based_boundary_id_offset = 0)
    {
      build_it(one_based_region_id,
               box_half_width,
               box_half_length,
               offset,
               one_based_boundary_id_offset);
    }


  private:
    /// Build the thing
    void build_it(const unsigned& one_based_region_id,
                  const double& box_half_width,
                  const double& box_half_length,
                  const Vector<double>& offset,
                  const unsigned& one_based_boundary_id_offset)
    {
      // Make vertices
      unsigned n_vertex = 8;
      Vertex_pt.resize(n_vertex);
      Vector<double> box_point(3);

      box_point[0] = offset[0] - box_half_width;
      box_point[1] = offset[1] - box_half_width;
      box_point[2] = offset[2] - box_half_length;
      Vertex_pt[0] = new TetMeshVertex(box_point);

      box_point[0] = offset[0] - box_half_width;
      box_point[1] = offset[1] + box_half_width;
      box_point[2] = offset[2] - box_half_length;
      Vertex_pt[1] = new TetMeshVertex(box_point);

      box_point[0] = offset[0] - box_half_width;
      box_point[1] = offset[1] + box_half_width;
      box_point[2] = offset[2] + box_half_length;
      Vertex_pt[2] = new TetMeshVertex(box_point);

      box_point[0] = offset[0] - box_half_width;
      box_point[1] = offset[1] - box_half_width;
      box_point[2] = offset[2] + box_half_length;
      Vertex_pt[3] = new TetMeshVertex(box_point);

      box_point[0] = offset[0] + box_half_width;
      box_point[1] = offset[1] - box_half_width;
      box_point[2] = offset[2] - box_half_length;
      Vertex_pt[4] = new TetMeshVertex(box_point);

      box_point[0] = offset[0] + box_half_width;
      box_point[1] = offset[1] + box_half_width;
      box_point[2] = offset[2] - box_half_length;
      Vertex_pt[5] = new TetMeshVertex(box_point);

      box_point[0] = offset[0] + box_half_width;
      box_point[1] = offset[1] + box_half_width;
      box_point[2] = offset[2] + box_half_length;
      Vertex_pt[6] = new TetMeshVertex(box_point);

      box_point[0] = offset[0] + box_half_width;
      box_point[1] = offset[1] - box_half_width;
      box_point[2] = offset[2] + box_half_length;
      Vertex_pt[7] = new TetMeshVertex(box_point);


      // Make facets
      unsigned n_facet = 6;
      Facet_pt.resize(n_facet);

      unsigned n_vertex_on_facet = 4;
      Facet_pt[0] = new TetMeshFacet(n_vertex_on_facet);
      unsigned one_based_boundary_id = 1 + one_based_boundary_id_offset;
      Facet_pt[0]->set_one_based_boundary_id(one_based_boundary_id);
      Facet_pt[0]->set_vertex_pt(0, Vertex_pt[0]);
      Facet_pt[0]->set_vertex_pt(1, Vertex_pt[4]);
      Facet_pt[0]->set_vertex_pt(2, Vertex_pt[7]);
      Facet_pt[0]->set_vertex_pt(3, Vertex_pt[3]);

      Facet_pt[1] = new TetMeshFacet(n_vertex_on_facet);
      one_based_boundary_id = 2 + one_based_boundary_id_offset;
      Facet_pt[1]->set_one_based_boundary_id(one_based_boundary_id);
      Facet_pt[1]->set_vertex_pt(0, Vertex_pt[4]);
      Facet_pt[1]->set_vertex_pt(1, Vertex_pt[5]);
      Facet_pt[1]->set_vertex_pt(2, Vertex_pt[6]);
      Facet_pt[1]->set_vertex_pt(3, Vertex_pt[7]);

      // top
      Facet_pt[2] = new TetMeshFacet(n_vertex_on_facet);
      one_based_boundary_id = 3 + one_based_boundary_id_offset;
      Facet_pt[2]->set_one_based_boundary_id(one_based_boundary_id);
      Facet_pt[2]->set_vertex_pt(0, Vertex_pt[3]);
      Facet_pt[2]->set_vertex_pt(1, Vertex_pt[7]);
      Facet_pt[2]->set_vertex_pt(2, Vertex_pt[6]);
      Facet_pt[2]->set_vertex_pt(3, Vertex_pt[2]);

      Facet_pt[3] = new TetMeshFacet(n_vertex_on_facet);
      one_based_boundary_id = 4 + one_based_boundary_id_offset;
      Facet_pt[3]->set_one_based_boundary_id(one_based_boundary_id);
      Facet_pt[3]->set_vertex_pt(0, Vertex_pt[6]);
      Facet_pt[3]->set_vertex_pt(1, Vertex_pt[5]);
      Facet_pt[3]->set_vertex_pt(2, Vertex_pt[1]);
      Facet_pt[3]->set_vertex_pt(3, Vertex_pt[2]);

      Facet_pt[4] = new TetMeshFacet(n_vertex_on_facet);
      one_based_boundary_id = 5 + one_based_boundary_id_offset;
      Facet_pt[4]->set_one_based_boundary_id(one_based_boundary_id);
      Facet_pt[4]->set_vertex_pt(0, Vertex_pt[3]);
      Facet_pt[4]->set_vertex_pt(1, Vertex_pt[2]);
      Facet_pt[4]->set_vertex_pt(2, Vertex_pt[1]);
      Facet_pt[4]->set_vertex_pt(3, Vertex_pt[0]);

      // bottom
      Facet_pt[5] = new TetMeshFacet(n_vertex_on_facet);
      one_based_boundary_id = 6 + one_based_boundary_id_offset;
      Facet_pt[5]->set_one_based_boundary_id(one_based_boundary_id);
      Facet_pt[5]->set_vertex_pt(0, Vertex_pt[5]);
      Facet_pt[5]->set_vertex_pt(1, Vertex_pt[1]);
      Facet_pt[5]->set_vertex_pt(2, Vertex_pt[0]);
      Facet_pt[5]->set_vertex_pt(3, Vertex_pt[4]);


      // Specify the adjacent region that is bounded by these facets
      for (unsigned f = 0; f < 6; f++)
      {
        Facet_pt[f]->set_one_based_adjacent_region_id(one_based_region_id);
      }

      // Is it a genuine region?
      if (one_based_region_id > 0)
      {
        set_region_for_tetgen(one_based_region_id - 1, offset);
      }
      // Otherwise declare it to be a hole
      else
      {
        set_hole_for_tetgen(offset);
        enable_faceted_volume_represents_hole_for_gmsh();
      }
    }
  };

  //=============================================================================
  /// RectangularTetMeshFacetedSurface from (+/- half x width) x +/- half y
  /// length
  //=============================================================================
  class RectangularTetMeshFacetedSurface : public virtual TetMeshFacetedSurface
  {
  public:
    /// Constructor: Specify dimension, offset vector and
    /// (one-based!) boundary_id.
    RectangularTetMeshFacetedSurface(const double& half_x_width,
                                     const double& half_y_length,
                                     const Vector<double>& offset,
                                     const unsigned& one_based_boundary_id)
    {
      build_it(half_x_width, half_y_length, offset, one_based_boundary_id);
    }


  private:
    /// Build the thing
    void build_it(const double& half_x_width,
                  const double& half_y_length,
                  const Vector<double>& offset,
                  const unsigned& one_based_boundary_id)
    {
      // Make vertices
      unsigned n_vertex = 4;
      Vertex_pt.resize(n_vertex);
      Vector<double> box_point(3);

      box_point[0] = offset[0] - half_x_width;
      box_point[1] = offset[1] - half_y_length;
      box_point[2] = offset[2];
      Vertex_pt[0] = new TetMeshVertex(box_point);

      box_point[0] = offset[0] + half_x_width;
      box_point[1] = offset[1] - half_y_length;
      box_point[2] = offset[2];
      Vertex_pt[1] = new TetMeshVertex(box_point);

      box_point[0] = offset[0] + half_x_width;
      box_point[1] = offset[1] + half_y_length;
      box_point[2] = offset[2];
      Vertex_pt[2] = new TetMeshVertex(box_point);

      box_point[0] = offset[0] - half_x_width;
      box_point[1] = offset[1] + half_y_length;
      box_point[2] = offset[2];
      Vertex_pt[3] = new TetMeshVertex(box_point);


      // Make facets
      unsigned n_facet = 1;
      Facet_pt.resize(n_facet);

      unsigned n_vertex_on_facet = 4;
      Facet_pt[0] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[0]->set_one_based_boundary_id(one_based_boundary_id);
      Facet_pt[0]->set_vertex_pt(0, Vertex_pt[0]);
      Facet_pt[0]->set_vertex_pt(1, Vertex_pt[1]);
      Facet_pt[0]->set_vertex_pt(2, Vertex_pt[2]);
      Facet_pt[0]->set_vertex_pt(3, Vertex_pt[3]);
    }
  };


  //================================================================
  /// TetMeshFacetedSurface that defines disk
  //================================================================
  class DiskTetMeshFacetedSurface : public virtual TetMeshFacetedSurface
  {
  public:
    /// Constructor: Pass pointer to GeomObject (with boundaries, and
    /// parametrised by coordinates without coordinate singularities, i.e. not
    /// polars, say) that defines the shape of the disk. Other args
    /// specify half the number of segments on perimeter of disk
    /// and first one-based boundary ID to be used to enumerate the boundaries
    /// on the disk. Returns last one-based boundary id used to enumerate the
    /// disk
    DiskTetMeshFacetedSurface(
      DiskLikeGeomObjectWithBoundaries*
        disk_parametrised_by_nonsingular_coordinates_pt,
      const unsigned& half_nsegment,
      const unsigned& first_one_based_boundary_id_for_disk,
      unsigned& last_one_based_boundary_id_for_disk)

    {
      Geom_object_with_boundaries_pt =
        disk_parametrised_by_nonsingular_coordinates_pt;

      // Provide storage for pointers to the two parts of the curvilinear
      // boundary
      Vector<TriangleMeshCurveSection*> outer_curvilinear_boundary_pt(2);

      // First bit
      GeomObject* outer_boundary_ellipse0_pt =
        disk_parametrised_by_nonsingular_coordinates_pt
          ->boundary_parametrising_geom_object_pt(0);
      double zeta_start =
        disk_parametrised_by_nonsingular_coordinates_pt->zeta_boundary_start(0);
      double zeta_end =
        disk_parametrised_by_nonsingular_coordinates_pt->zeta_boundary_end(0);
      unsigned nsegment = half_nsegment;
      unsigned boundary_id = 0;
      outer_curvilinear_boundary_pt[0] =
        new TriangleMeshCurviLine(outer_boundary_ellipse0_pt,
                                  zeta_start,
                                  zeta_end,
                                  nsegment,
                                  boundary_id);

      // Second bit
      GeomObject* outer_boundary_ellipse1_pt =
        disk_parametrised_by_nonsingular_coordinates_pt
          ->boundary_parametrising_geom_object_pt(1);
      zeta_start =
        disk_parametrised_by_nonsingular_coordinates_pt->zeta_boundary_start(1);
      zeta_end =
        disk_parametrised_by_nonsingular_coordinates_pt->zeta_boundary_end(1);
      nsegment = half_nsegment;
      boundary_id = 1;
      outer_curvilinear_boundary_pt[1] =
        new TriangleMeshCurviLine(outer_boundary_ellipse1_pt,
                                  zeta_start,
                                  zeta_end,
                                  nsegment,
                                  boundary_id);

      // Combine to curvilinear boundary and define the
      // outer boundary
      TriangleMeshClosedCurve* closed_curve_pt =
        new TriangleMeshClosedCurve(outer_curvilinear_boundary_pt);

      // Use the TriangleMeshParameters object for helping on the manage of the
      // TriangleMesh parameters
      TriangleMeshParameters triangle_mesh_parameters(closed_curve_pt);

      // Do we have an annnular internal boundary?
      if (disk_parametrised_by_nonsingular_coordinates_pt->nboundary() == 4)
      {
        // Provide storage for pointers to the two parts of the curvilinear
        // boundary
        Vector<TriangleMeshCurveSection*> inner_curvilinear_boundary_pt(2);

        // First bit
        GeomObject* inner_boundary_ellipse0_pt =
          disk_parametrised_by_nonsingular_coordinates_pt
            ->boundary_parametrising_geom_object_pt(2);
        double zeta_start =
          disk_parametrised_by_nonsingular_coordinates_pt->zeta_boundary_start(
            2);
        double zeta_end =
          disk_parametrised_by_nonsingular_coordinates_pt->zeta_boundary_end(2);
        unsigned nsegment = half_nsegment;
        unsigned boundary_id = 2;
        inner_curvilinear_boundary_pt[0] =
          new TriangleMeshCurviLine(inner_boundary_ellipse0_pt,
                                    zeta_start,
                                    zeta_end,
                                    nsegment,
                                    boundary_id);

        // Second bit
        GeomObject* inner_boundary_ellipse1_pt =
          disk_parametrised_by_nonsingular_coordinates_pt
            ->boundary_parametrising_geom_object_pt(3);
        zeta_start =
          disk_parametrised_by_nonsingular_coordinates_pt->zeta_boundary_start(
            3);
        zeta_end =
          disk_parametrised_by_nonsingular_coordinates_pt->zeta_boundary_end(3);
        nsegment = half_nsegment;
        boundary_id = 3;
        inner_curvilinear_boundary_pt[1] =
          new TriangleMeshCurviLine(inner_boundary_ellipse1_pt,
                                    zeta_start,
                                    zeta_end,
                                    nsegment,
                                    boundary_id);

        // Combine to curvilinear boundary and define the
        // inner boundary
        TriangleMeshClosedCurve* inner_curve_pt =
          new TriangleMeshClosedCurve(inner_curvilinear_boundary_pt);


        Vector<TriangleMeshClosedCurve*> inner_boundaries_pt(1);
        inner_boundaries_pt[0] = inner_curve_pt;

        // Specify the internal closed boundaries
        triangle_mesh_parameters.internal_closed_curve_pt() =
          inner_boundaries_pt;


        // Let's get a map of region coordinates (if any!)
        std::map<unsigned, Vector<double>> zeta_in_region_map =
          disk_parametrised_by_nonsingular_coordinates_pt->zeta_in_region();
        for (std::map<unsigned, Vector<double>>::iterator it =
               zeta_in_region_map.begin();
             it != zeta_in_region_map.end();
             it++)
        {
          // Pass information about the defined regions
          unsigned r = (*it).first;
          Vector<double> region_coords = (*it).second;
          triangle_mesh_parameters.add_region_coordinates(r, region_coords);
        }
      }

      // Specify the element area so that it matches the boundary discretisation
      double uniform_element_area =
        sqrt(3.0) / 4.0 *
        pow(2.0 * MathematicalConstants::Pi / (2.0 * double(half_nsegment)), 2);
      triangle_mesh_parameters.element_area() = uniform_element_area;

      // Create the mesh
      Tri_mesh_pt =
        new TriangleMesh<TPoissonElement<2, 2>>(triangle_mesh_parameters);

      // Loop over all boundary elements and rotate their nodes so that
      // the first two nodes are on the boundary
      std::map<FiniteElement*, bool> is_on_boundary;
      for (unsigned b = 0; b < 2; b++)
      {
        unsigned nel = Tri_mesh_pt->nboundary_element(b);
        for (unsigned e = 0; e < nel; e++)
        {
          FiniteElement* el_pt = Tri_mesh_pt->boundary_element_pt(b, e);
          unsigned count = 0;
          for (unsigned j = 0; j < 3; j++)
          {
            Node* nod_pt = el_pt->node_pt(j);
            if (nod_pt->is_on_boundary()) count++;
          }
          if (count == 2)
          {
            is_on_boundary[el_pt] = true;
            if ((el_pt->node_pt(0)->is_on_boundary()) &&
                (el_pt->node_pt(1)->is_on_boundary()))
            {
              // fine
            }
            else
            {
              // Reorder nodes so that the first two nodes are on the boundary
              Node* nod0_pt = el_pt->node_pt(0);
              Node* nod1_pt = el_pt->node_pt(1);
              Node* nod2_pt = el_pt->node_pt(2);
              if (!el_pt->node_pt(0)->is_on_boundary())
              {
                el_pt->node_pt(0) = nod1_pt;
                el_pt->node_pt(1) = nod2_pt;
                el_pt->node_pt(2) = nod0_pt;
              }
              else if (!el_pt->node_pt(1)->is_on_boundary())
              {
                el_pt->node_pt(0) = nod2_pt;
                el_pt->node_pt(1) = nod0_pt;
                el_pt->node_pt(2) = nod1_pt;
              }
              else
              {
                std::ostringstream error_message;
                error_message << "This doesn't make sense!";
                throw OomphLibError(error_message.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
            }
          }
          else
          {
            std::ostringstream error_message;
            error_message
              << "Boundary simplex element doesn't have two nodes on boundary!";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
      }


      // Now loop over all nodes and turn them into vertices
      unsigned nnod = Tri_mesh_pt->nnode();
      Vertex_pt.resize(nnod);
      for (unsigned j = 0; j < nnod; j++)
      {
        Node* nod_pt = Tri_mesh_pt->node_pt(j);
        Vector<double> sheet_point(3, 0.0);
        Vector<double> zeta(nod_pt->position());
        Geom_object_with_boundaries_pt->position(zeta, sheet_point);
        Vertex_pt[j] = new TetMeshVertex(sheet_point);
        Vertex_pt[j]->set_zeta_in_geom_object(zeta);
        Equivalent_vertex_pt[nod_pt] = Vertex_pt[j];
      }


      // Quick count to see how many elements/facets are on the boundary
      Nelement_on_disk_boundary = 0;
      unsigned nel = Tri_mesh_pt->nelement();
      for (unsigned e = 0; e < nel; e++)
      {
        FiniteElement* el_pt = Tri_mesh_pt->finite_element_pt(e);
        if (is_on_boundary[el_pt])
        {
          Nelement_on_disk_boundary++;
        }
      }

      // Angles of right/left node on boundary (for snapping)
      Left_boundary_coordinate.resize(nel, 0.0);
      Right_boundary_coordinate.resize(nel, 0.0);
      Boundary_id.resize(nel, 0.0);

      // Now loop over all elements
      Nfacet = nel;
      Facet_pt.resize(Nfacet);
      Facet_is_on_boundary.resize(Nfacet, false);
      unsigned n_vertex_on_facet = 3;
      for (unsigned e = 0; e < nel; e++)
      {
        FiniteElement* el_pt = Tri_mesh_pt->finite_element_pt(e);
        if (is_on_boundary[el_pt])
        {
          Node* left_node_pt = el_pt->node_pt(0);
          Node* right_node_pt = el_pt->node_pt(1);
          Vector<double> boundary_zeta(1);
          if (left_node_pt->is_on_boundary(0) &&
              right_node_pt->is_on_boundary(0))
          {
            unsigned b = 0;
            left_node_pt->get_coordinates_on_boundary(b, boundary_zeta);
            Left_boundary_coordinate[e] = boundary_zeta[0];
            right_node_pt->get_coordinates_on_boundary(b, boundary_zeta);
            Right_boundary_coordinate[e] = boundary_zeta[0];
            Boundary_id[e] = b;
          }
          else if (left_node_pt->is_on_boundary(1) &&
                   right_node_pt->is_on_boundary(1))
          {
            unsigned b = 1;
            left_node_pt->get_coordinates_on_boundary(b, boundary_zeta);
            Left_boundary_coordinate[e] = boundary_zeta[0];
            right_node_pt->get_coordinates_on_boundary(b, boundary_zeta);
            Right_boundary_coordinate[e] = boundary_zeta[0];
            Boundary_id[e] = b;
          }
          else
          {
            std::ostringstream error_message;
            error_message << "never get here!";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
          Facet_is_on_boundary[e] = true;
        }
        Facet_pt[e] = new TetMeshFacet(n_vertex_on_facet);
        unsigned one_based_boundary_id =
          first_one_based_boundary_id_for_disk + e;
        Facet_pt[e]->set_one_based_boundary_id(one_based_boundary_id);
        last_one_based_boundary_id_for_disk = one_based_boundary_id;
        for (unsigned j = 0; j < 3; j++)
        {
          Facet_pt[e]->set_vertex_pt(j,
                                     Equivalent_vertex_pt[el_pt->node_pt(j)]);
        }
      }
    }


    /// Destructor
    ~DiskTetMeshFacetedSurface()
    {
      delete Tri_mesh_pt;
      Tri_mesh_pt = 0;
    }

    /// Function that specifies the variation of the
    /// zeta coordinates in the GeomObject along the
    /// boundary connecting vertices 0 and 1 in the
    /// facet. NOTE: input zeta_boundary ranges between 0 and 1;
    /// gets mapped to actual boundary coordinate inside!
    void boundary_zeta01(const unsigned& facet_id,
                         const double& zeta_boundary,
                         Vector<double>& zeta)
    {
      if (Facet_is_on_boundary[facet_id])
      {
        double phi0 = Left_boundary_coordinate[facet_id];
        double phi1 = Right_boundary_coordinate[facet_id];
        double phi = phi0 + (phi1 - phi0) * zeta_boundary;

        unsigned b = Boundary_id[facet_id];
        double zeta_bound = phi;
        Geom_object_with_boundaries_pt->zeta_on_boundary(b, zeta_bound, zeta);
      }
      else
      {
        Vector<double> zeta0 =
          Facet_pt[facet_id]->vertex_pt(0)->zeta_in_geom_object();
        Vector<double> zeta1 =
          Facet_pt[facet_id]->vertex_pt(1)->zeta_in_geom_object();
        zeta[0] = zeta0[0] + (zeta1[0] - zeta0[0]) * zeta_boundary;
        zeta[1] = zeta0[1] + (zeta1[1] - zeta0[1]) * zeta_boundary;
      }
    }


  protected:
    /// Number of facets
    unsigned Nfacet;

    /// Is facet on boundary?
    std::vector<bool> Facet_is_on_boundary;

    /// Left boundary coordinate of i-th facet
    Vector<double> Left_boundary_coordinate;

    /// Right boundary coordinate of i-th facet
    Vector<double> Right_boundary_coordinate;

    /// ID of boundary the i-th facet is located on
    Vector<unsigned> Boundary_id;

    /// Number of elements on disk boundary
    unsigned Nelement_on_disk_boundary;


    /// Mesh used to facet-ise (discretise) disk
    TriangleMesh<TPoissonElement<2, 2>>* Tri_mesh_pt;


    /// Mapping between nodes and vertices
    std::map<Node*, TetMeshVertex*> Equivalent_vertex_pt;
  };


  //================================================================
  /// TetMeshFacetedSurface that defines disk with torus around edge
  //================================================================
  class DiskWithTorusAroundEdgeTetMeshFacetedSurface
    : public virtual DiskTetMeshFacetedSurface,
      public virtual TetMeshFacetedClosedSurface
  {
  public:
    /// Constructor: Pass pointer to GeomObject (with boundaries, and
    /// parametrised by coordinates without coordinate singularities, i.e. not
    /// polars, say) that defines the shape of the disk. Other args
    /// specify:
    /// - half the number of segments on perimeter of disk
    /// - radius of torus
    /// - number of vertices around the perimeter of the torus
    /// - first one-based boundary ID to be used to enumerate the boundaries on
    /// the
    ///   disk
    /// - one based region ID for volume contained inside torus.
    /// .
    /// Computes:
    /// - last one-based boundary id used to enumerate boundaries on disk
    /// - first one-based boundary id used to enumerate boundaries on torus
    /// - last one-based boundary id used to enumerate boundaries on torus
    /// - vector containing one-based boundary ids of boundaries on disk
    ///   enclosed within the torus
    /// - vector containing one-based boundary ids of boundaries on disk
    ///   not enclosed within the torus
    /// .
    DiskWithTorusAroundEdgeTetMeshFacetedSurface(
      DiskLikeGeomObjectWithBoundaries*
        disk_parametrised_by_nonsingular_coordinates_pt,
      const unsigned& half_nsegment,
      const double& r_torus,
      const unsigned& nvertex_torus,
      const unsigned& first_one_based_boundary_id_for_disk,
      const unsigned& one_based_torus_region_id,
      unsigned& last_one_based_boundary_id_for_disk,
      unsigned& first_one_based_boundary_id_for_torus,
      unsigned& last_one_based_boundary_id_for_torus,
      Vector<unsigned>& one_based_boundary_id_for_disk_within_torus,
      Vector<unsigned>& one_based_boundary_id_for_disk_outside_torus)
      : DiskTetMeshFacetedSurface(
          disk_parametrised_by_nonsingular_coordinates_pt,
          half_nsegment,
          first_one_based_boundary_id_for_disk,
          last_one_based_boundary_id_for_disk)
    {
      // Wipe
      one_based_boundary_id_for_disk_within_torus.clear();
      one_based_boundary_id_for_disk_outside_torus.clear();


      // Is element in torus region or not?
      std::map<FiniteElement*, bool> is_in_torus_region;
      {
        unsigned r = 1;
        /* ofstream some_file; */
        /* some_file.open("tri_mesh_in_region1.dat");  */
        unsigned nel = Tri_mesh_pt->nregion_element(r);
        for (unsigned e = 0; e < nel; e++)
        {
          FiniteElement* el_pt = Tri_mesh_pt->region_element_pt(r, e);
          is_in_torus_region[el_pt] = true;
          // unsigned npts=3;
          // el_pt->output(some_file,npts);
        }
        // some_file.close();
      }


      // Now loop over all elements
      unsigned nel = Tri_mesh_pt->nelement();
      for (unsigned e = 0; e < nel; e++)
      {
        FiniteElement* el_pt = Tri_mesh_pt->finite_element_pt(e);
        unsigned one_based_boundary_id = Facet_pt[e]->one_based_boundary_id();
        if (is_in_torus_region[el_pt])
        {
          Facet_pt[e]->set_one_based_region_that_facet_is_embedded_in(
            one_based_torus_region_id);
          one_based_boundary_id_for_disk_within_torus.push_back(
            one_based_boundary_id);
        }
        else
        {
          one_based_boundary_id_for_disk_outside_torus.push_back(
            one_based_boundary_id);
        }
      }


      // Now add torus:
      //---------------
      unsigned one_based_boundary_id = last_one_based_boundary_id_for_disk + 1;
      first_one_based_boundary_id_for_torus = one_based_boundary_id;

      // Find sorted sequence of vertices going around the
      // inner boundary (innermost part of torus)
      Vector<TetMeshVertex*> sorted_vertex_pt;
      Vector<std::pair<TetMeshVertex*, double>> vertex_pt_and_boundary_coord;
      std::map<Node*, bool> done;
      Vector<double> boundary_zeta(1);
      for (unsigned b = 2; b <= 3; b++)
      {
        unsigned nnod = Tri_mesh_pt->nboundary_node(b);
        for (unsigned j = 0; j < nnod; j++)
        {
          Node* nod_pt = Tri_mesh_pt->boundary_node_pt(b, j);
          if (!done[nod_pt])
          {
            nod_pt->get_coordinates_on_boundary(b, boundary_zeta);
            vertex_pt_and_boundary_coord.push_back(
              std::make_pair(Equivalent_vertex_pt[nod_pt], boundary_zeta[0]));
            done[nod_pt] = true;
          }
        }
      }


      oomph_info << "Number of vertices on inner ring: "
                 << vertex_pt_and_boundary_coord.size() << " "
                 << " sum: "
                 << Tri_mesh_pt->nboundary_node(2) +
                      Tri_mesh_pt->nboundary_node(3)
                 << std::endl;

      // Identify point in torus (for tetgen)
      Vector<double> point_in_torus(3);
      unsigned b = 0;
      unsigned j = 0;
      Node* nod_pt = Tri_mesh_pt->boundary_node_pt(b, j);
      point_in_torus[0] = Equivalent_vertex_pt[nod_pt]->x(0);
      point_in_torus[1] = Equivalent_vertex_pt[nod_pt]->x(1);
      point_in_torus[2] = Equivalent_vertex_pt[nod_pt]->x(2) + 0.5 * r_torus;

      oomph_info << "point in torus: " << point_in_torus[0] << " "
                 << point_in_torus[1] << " " << point_in_torus[2] << " "
                 << std::endl;

      // Sort based on boundary coordinate
      std::sort(vertex_pt_and_boundary_coord.begin(),
                vertex_pt_and_boundary_coord.end(),
                HelperForSortingVectorOfPairsOfVertexPtAndDouble::
                  less_than_based_on_double);

      // Loop over annular rings
      unsigned n_vertices_on_annular_ring = nvertex_torus;
      unsigned n_new_vertices_on_annular_ring = n_vertices_on_annular_ring - 2;


      // Storage for vertices and facets
      unsigned nv = Vertex_pt.size();
      unsigned jv = nv;
      nv +=
        n_new_vertices_on_annular_ring * vertex_pt_and_boundary_coord.size();
      Vertex_pt.resize(nv);

      unsigned jf = Nfacet;
      Nfacet += (n_new_vertices_on_annular_ring + 1) *
                (vertex_pt_and_boundary_coord.size());
      Facet_pt.resize(Nfacet);


      unsigned new_vertex_count = 0;
      unsigned new_facet_count = 0;
      unsigned n = vertex_pt_and_boundary_coord.size();
      Vector<Vector<TetMeshVertex*>> annulus_vertex_pt(n);
      for (unsigned j = 0; j < n; j++)
      {
        // Wrap around for simplicity
        annulus_vertex_pt[j].resize(n_vertices_on_annular_ring);
        annulus_vertex_pt[j][0] = vertex_pt_and_boundary_coord[j].first;
        annulus_vertex_pt[j][n_vertices_on_annular_ring - 1] =
          vertex_pt_and_boundary_coord[j].first;

        // Position of vertex
        Vector<double> vertex_position(3);
        vertex_position[0] = vertex_pt_and_boundary_coord[j].first->x(0);
        vertex_position[1] = vertex_pt_and_boundary_coord[j].first->x(1);
        vertex_position[2] = vertex_pt_and_boundary_coord[j].first->x(2);

        // Polar angle
        double theta = vertex_pt_and_boundary_coord[j].second;

        // Which boundary is this point located on?
        unsigned b = 0;
        if (theta > Geom_object_with_boundaries_pt->zeta_boundary_end(0))
        {
          b = 1;
        }
        double zeta_bound = theta;

        // Get coordinate on edge and normal and tangent vectors
        Vector<double> x(3);
        Vector<double> r_edge(3);
        Vector<double> normal(3);
        Vector<double> tangent(3);
        Vector<double> normal_normal(3);
        Geom_object_with_boundaries_pt->boundary_triad(
          b, zeta_bound, r_edge, tangent, normal, normal_normal);

        // What rho do we need to create a circle, centred at the edge
        // going through current vertex?
        Vector<double> distance_to_edge(3);
        distance_to_edge[0] = r_edge[0] - vertex_position[0];
        distance_to_edge[1] = r_edge[1] - vertex_position[1];
        distance_to_edge[2] = r_edge[2] - vertex_position[2];
        double rho = sqrt(distance_to_edge[0] * distance_to_edge[0] +
                          distance_to_edge[1] * distance_to_edge[1] +
                          distance_to_edge[2] * distance_to_edge[2]);

        // What is the starting angle?
        double cos_phi0 =
          (distance_to_edge[0] * normal[0] + distance_to_edge[1] * normal[1] +
           distance_to_edge[2] * normal[2]) /
          rho;
        if (cos_phi0 > 1.0) cos_phi0 = 1.0;
        if (cos_phi0 < -1.0) cos_phi0 = -1.0;

        double phi_0 = acos(cos_phi0);

        for (unsigned i = 1; i < n_vertices_on_annular_ring - 1; i++)
        {
          double phi = phi_0 +
                       double(i) / double(n_vertices_on_annular_ring - 1) *
                         2.0 * MathematicalConstants::Pi -
                       MathematicalConstants::Pi;
          x[0] = r_edge[0] + rho * cos(phi) * normal[0] +
                 rho * sin(phi) * normal_normal[0];
          x[1] = r_edge[1] + rho * cos(phi) * normal[1] +
                 rho * sin(phi) * normal_normal[1];
          x[2] = r_edge[2] + rho * cos(phi) * normal[2] +
                 rho * sin(phi) * normal_normal[2];

          TetMeshVertex* new_vertex_pt = new TetMeshVertex(x);
          annulus_vertex_pt[j][i] = new_vertex_pt;
          Vertex_pt[jv] = new_vertex_pt;
          jv++;
          new_vertex_count++;
        }
      }

      for (unsigned i = 0; i < n - 1; i++)
      {
        unsigned m = annulus_vertex_pt[i].size();
        for (unsigned j = 0; j < m - 1; j++)
        {
          unsigned n_vertex_on_facet = 4;
          TetMeshFacet* new_facet_pt = new TetMeshFacet(n_vertex_on_facet);
          new_facet_count++;
          new_facet_pt->set_vertex_pt(0, annulus_vertex_pt[i][j]);
          new_facet_pt->set_vertex_pt(1, annulus_vertex_pt[i][j + 1]);
          new_facet_pt->set_vertex_pt(2, annulus_vertex_pt[i + 1][j + 1]);
          new_facet_pt->set_vertex_pt(3, annulus_vertex_pt[i + 1][j]);
          Facet_pt[jf] = new_facet_pt;
          unsigned region_id = one_based_torus_region_id;
          Facet_pt[jf]->set_one_based_adjacent_region_id(region_id);
          Facet_pt[jf]->set_one_based_boundary_id(one_based_boundary_id);
          one_based_boundary_id++;
          jf++;
        }
      }

      unsigned m = annulus_vertex_pt[n - 1].size();
      for (unsigned j = 0; j < m - 1; j++)
      {
        unsigned n_vertex_on_facet = 4;
        TetMeshFacet* new_facet_pt = new TetMeshFacet(n_vertex_on_facet);
        new_facet_count++;
        new_facet_pt->set_vertex_pt(0, annulus_vertex_pt[n - 1][j]);
        new_facet_pt->set_vertex_pt(1, annulus_vertex_pt[n - 1][j + 1]);
        new_facet_pt->set_vertex_pt(2, annulus_vertex_pt[0][j + 1]);
        new_facet_pt->set_vertex_pt(3, annulus_vertex_pt[0][j]);
        Facet_pt[jf] = new_facet_pt;
        unsigned region_id = one_based_torus_region_id;
        Facet_pt[jf]->set_one_based_adjacent_region_id(region_id);
        Facet_pt[jf]->set_one_based_boundary_id(one_based_boundary_id);
        last_one_based_boundary_id_for_torus = one_based_boundary_id;
        one_based_boundary_id++;
        jf++;
      }

      oomph_info << "First/last one-based disk boundary ID: "
                 << first_one_based_boundary_id_for_disk << " "
                 << last_one_based_boundary_id_for_disk << std::endl;


      oomph_info << "First/last one-based torus boundary ID: "
                 << first_one_based_boundary_id_for_torus << " "
                 << last_one_based_boundary_id_for_torus << std::endl;


      // Is it a genuine region?
      if (one_based_torus_region_id > 0)
      {
        set_region_for_tetgen(one_based_torus_region_id - 1, point_in_torus);
      }
      else
      {
        std::ostringstream error_message;
        error_message << "one_based_torus_region_id must be strictly positive";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
  };


  //================================================================
  /// TetMeshFacetedClosedSurface that defines disk with two layers
  /// on top and bottom
  //================================================================
  class DiskWithTwoLayersTetMeshFacetedSurface
    : public virtual DiskTetMeshFacetedSurface,
      public virtual TetMeshFacetedClosedSurface
  {
  public:
    /// Constructor: Pass pointer to GeomObject (with boundaries, and
    /// parametrised by coordinates without coordinate singularities, i.e. not
    /// polars, say) that defines the shape of the disk. Other args
    /// specify:
    /// - half the number of segments on perimeter of disk
    /// - first one-based boundary ID to be used to enumerate the boundaries on
    /// the
    ///   disk
    /// - one based region ID for region above the disk
    /// - one based region ID for region below the disk
    /// .
    /// Computes:
    /// - last one-based boundary id used to enumerate boundaries on disk
    /// - last one-based boundary id used to enumerate any boundaries
    /// .
    DiskWithTwoLayersTetMeshFacetedSurface(
      DiskLikeGeomObjectWithBoundaries*
        disk_parametrised_by_nonsingular_coordinates_pt,
      const unsigned& half_nsegment,
      const unsigned& first_one_based_boundary_id_for_disk,
      const unsigned& one_based_region_id_above_disk,
      const unsigned& one_based_region_id_below_disk,
      unsigned& last_one_based_boundary_id_for_disk,
      unsigned& last_one_based_boundary_id)
      : DiskTetMeshFacetedSurface(
          disk_parametrised_by_nonsingular_coordinates_pt,
          half_nsegment,
          first_one_based_boundary_id_for_disk,
          last_one_based_boundary_id_for_disk)
    {
      // Currently: Number of facets = number of elements on disk.
      // Set bounding region numbers for existing facets:
      unsigned nel = Nfacet;
      for (unsigned e = 0; e < Nfacet; e++)
      {
        Facet_pt[e]->set_one_based_adjacent_region_id(
          one_based_region_id_above_disk);
        Facet_pt[e]->set_one_based_adjacent_region_id(
          one_based_region_id_below_disk);
      }

      // Additional facets: one square facet moving up/down from each
      // side and two lids.
      Nfacet += 2 * Nelement_on_disk_boundary + 2;
      Facet_pt.resize(Nfacet);
      Facet_is_on_boundary.resize(Nfacet, false);

      // Find sorted sequence of vertices going around the
      // perimeter of the disk
      Vector<TetMeshVertex*> sorted_perimeter_vertex_pt;
      std::set<TetMeshFacet*> boundary_facet_pt;
      bool first = true;
      for (unsigned e = 0; e < nel; e++)
      {
        if (Facet_is_on_boundary[e])
        {
          boundary_facet_pt.insert(Facet_pt[e]);
          if (first)
          {
            TetMeshVertex* vertex0_pt = Facet_pt[e]->vertex_pt(0);
            sorted_perimeter_vertex_pt.push_back(vertex0_pt);
            first = false;
          }
        }
      }
      while (boundary_facet_pt.size() != 1)
      {
        for (std::set<TetMeshFacet*>::iterator it = boundary_facet_pt.begin();
             it != boundary_facet_pt.end();
             it++)
        {
          unsigned n_last = sorted_perimeter_vertex_pt.size();
          if (sorted_perimeter_vertex_pt[n_last - 1] == (*it)->vertex_pt(0))
          {
            sorted_perimeter_vertex_pt.push_back((*it)->vertex_pt(1));
            boundary_facet_pt.erase(*it);
            break;
          }
        }
      }

      // Find mean position
      double z_mean = 0.0;
      double z_min = DBL_MAX;
      double z_max = -DBL_MAX;
      unsigned n = sorted_perimeter_vertex_pt.size();
      for (unsigned j = 0; j < n; j++)
      {
        double z = sorted_perimeter_vertex_pt[j]->x(2);
        z_mean += z;
        if (z > z_max) z_max = z;
        if (z < z_min) z_min = z;
      }
      z_mean /= double(n);
      double layer_thickness = 0.6 * (z_max - z_min);

      // Loop over upper/lower layer
      unsigned facet_count = nel;
      int offset_sign = -1;
      unsigned one_based_boundary_id = last_one_based_boundary_id_for_disk + 1;
      unsigned n_vertex_on_disk = Vertex_pt.size();
      unsigned n_vertex_on_disk_perimeter = sorted_perimeter_vertex_pt.size();
      for (unsigned i_offset = 0; i_offset < 2; i_offset++)
      {
        unsigned one_based_region_id = one_based_region_id_below_disk;
        if (i_offset == 1) one_based_region_id = one_based_region_id_above_disk;

        for (unsigned j = 0; j < n_vertex_on_disk_perimeter; j++)
        {
          Vector<double> posn(3);
          posn[0] = sorted_perimeter_vertex_pt[j]->x(0);
          posn[1] = sorted_perimeter_vertex_pt[j]->x(1);
          posn[2] = z_mean + double(offset_sign) * layer_thickness;
          Vertex_pt.push_back(new TetMeshVertex(posn));
        }
        unsigned count = 0;
        for (unsigned j = 0; j < n_vertex_on_disk_perimeter - 1; j++)
        {
          Facet_pt[facet_count] = new TetMeshFacet(4);
          Facet_pt[facet_count]->set_one_based_boundary_id(
            one_based_boundary_id);
          Facet_pt[facet_count]->set_vertex_pt(0,
                                               sorted_perimeter_vertex_pt[j]);
          Facet_pt[facet_count]->set_vertex_pt(
            1, sorted_perimeter_vertex_pt[j + 1]);
          Facet_pt[facet_count]->set_vertex_pt(
            2,
            Vertex_pt[n_vertex_on_disk + i_offset * n_vertex_on_disk_perimeter +
                      j + 1]);
          Facet_pt[facet_count]->set_vertex_pt(
            3,
            Vertex_pt[n_vertex_on_disk + i_offset * n_vertex_on_disk_perimeter +
                      j]);
          Facet_pt[facet_count]->set_one_based_adjacent_region_id(
            one_based_region_id);
          count++;
          facet_count++;
          one_based_boundary_id++;
        }
        Facet_pt[facet_count] = new TetMeshFacet(4);
        Facet_pt[facet_count]->set_one_based_boundary_id(one_based_boundary_id);
        Facet_pt[facet_count]->set_vertex_pt(0,
                                             sorted_perimeter_vertex_pt[count]);
        Facet_pt[facet_count]->set_vertex_pt(1, sorted_perimeter_vertex_pt[0]);
        Facet_pt[facet_count]->set_vertex_pt(
          2,
          Vertex_pt[n_vertex_on_disk + i_offset * n_vertex_on_disk_perimeter +
                    0]);
        Facet_pt[facet_count]->set_vertex_pt(
          3,
          Vertex_pt[n_vertex_on_disk + i_offset * n_vertex_on_disk_perimeter +
                    count]);
        Facet_pt[facet_count]->set_one_based_adjacent_region_id(
          one_based_region_id);


        // Closing lid
        facet_count++;
        one_based_boundary_id++;
        Facet_pt[facet_count] = new TetMeshFacet(n_vertex_on_disk_perimeter);
        Facet_pt[facet_count]->set_one_based_boundary_id(one_based_boundary_id);
        for (unsigned j = 0; j < n_vertex_on_disk_perimeter; j++)
        {
          Facet_pt[facet_count]->set_vertex_pt(
            j,
            Vertex_pt[n_vertex_on_disk + i_offset * n_vertex_on_disk_perimeter +
                      j]);
        }
        Facet_pt[facet_count]->set_one_based_adjacent_region_id(
          one_based_region_id);

        // Last one used (gets assigned twice but last one sticks!)
        last_one_based_boundary_id = one_based_boundary_id;

        // Bump
        facet_count++;
        one_based_boundary_id++;

        /// Identify region for tetgen
        Vector<double> region_point(3, 0.0);
        region_point[2] = z_mean + 0.5 * double(offset_sign) * layer_thickness;
        set_region_for_tetgen(one_based_region_id - 1, region_point);

        offset_sign *= -1;
      }
    }
  };

  //=============================================================
  /// TetMeshFacetedSurface that defines inner boundary
  //=============================================================
  class SphericalTetMeshFacetedSurface
    : public virtual TetMeshFacetedClosedSurface
  {
  public:
    /// Constructor
    SphericalTetMeshFacetedSurface()
    {
      // Golden ratio
      const double phi = 0.5 * (1.0 + sqrt(5.0));

      // Make vertices
      unsigned n_vertex = 12;
      Vertex_pt.resize(n_vertex);

      // Set basic icosahedron points
      Vector<double> icosa_point(3);

      icosa_point[0] = 0.0;
      icosa_point[1] = 1.0;
      icosa_point[2] = phi;
      Vertex_pt[0] = new TetMeshVertex(icosa_point);

      icosa_point[0] = 0.0;
      icosa_point[1] = -1.0;
      icosa_point[2] = phi;
      Vertex_pt[1] = new TetMeshVertex(icosa_point);

      icosa_point[0] = 0.0;
      icosa_point[1] = 1.0;
      icosa_point[2] = -phi;
      Vertex_pt[2] = new TetMeshVertex(icosa_point);

      icosa_point[0] = 0.0;
      icosa_point[1] = -1.0;
      icosa_point[2] = -phi;
      Vertex_pt[3] = new TetMeshVertex(icosa_point);

      icosa_point[0] = 1.0;
      icosa_point[1] = phi;
      icosa_point[2] = 0.0;
      Vertex_pt[4] = new TetMeshVertex(icosa_point);

      icosa_point[0] = -1.0;
      icosa_point[1] = phi;
      icosa_point[2] = 0.0;
      Vertex_pt[5] = new TetMeshVertex(icosa_point);

      icosa_point[0] = 1.0;
      icosa_point[1] = -phi;
      icosa_point[2] = 0.0;
      Vertex_pt[6] = new TetMeshVertex(icosa_point);

      icosa_point[0] = -1.0;
      icosa_point[1] = -phi;
      icosa_point[2] = 0.0;
      Vertex_pt[7] = new TetMeshVertex(icosa_point);

      icosa_point[0] = phi;
      icosa_point[1] = 0.0;
      icosa_point[2] = 1.0;
      Vertex_pt[8] = new TetMeshVertex(icosa_point);

      icosa_point[0] = phi;
      icosa_point[1] = 0.0;
      icosa_point[2] = -1.0;
      Vertex_pt[9] = new TetMeshVertex(icosa_point);

      icosa_point[0] = -phi;
      icosa_point[1] = 0.0;
      icosa_point[2] = 1.0;
      Vertex_pt[10] = new TetMeshVertex(icosa_point);

      icosa_point[0] = -phi;
      icosa_point[1] = 0.0;
      icosa_point[2] = -1.0;
      Vertex_pt[11] = new TetMeshVertex(icosa_point);

      // Make facets
      unsigned n_facet = 20;
      Facet_pt.resize(n_facet);

      unsigned n_vertex_on_facet = 3;
      Facet_pt[0] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[0]->set_vertex_pt(0, Vertex_pt[0]);
      Facet_pt[0]->set_vertex_pt(1, Vertex_pt[1]);
      Facet_pt[0]->set_vertex_pt(2, Vertex_pt[8]);

      // icosa_facet[0][0] = 0;
      // icosa_facet[0][1] = 1;
      // icosa_facet[0][2] = 8;

      Facet_pt[1] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[1]->set_vertex_pt(0, Vertex_pt[0]);
      Facet_pt[1]->set_vertex_pt(1, Vertex_pt[10]);
      Facet_pt[1]->set_vertex_pt(2, Vertex_pt[1]);

      // icosa_facet[1].resize(3);
      // icosa_facet[1][0] = 0;
      // icosa_facet[1][1] = 10;
      // icosa_facet[1][2] = 1;

      Facet_pt[2] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[2]->set_vertex_pt(0, Vertex_pt[0]);
      Facet_pt[2]->set_vertex_pt(1, Vertex_pt[5]);
      Facet_pt[2]->set_vertex_pt(2, Vertex_pt[10]);

      // icosa_facet[2].resize(3);
      // icosa_facet[2][0] = 0;
      // icosa_facet[2][1] = 5;
      // icosa_facet[2][2] = 10;

      Facet_pt[3] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[3]->set_vertex_pt(0, Vertex_pt[0]);
      Facet_pt[3]->set_vertex_pt(1, Vertex_pt[4]);
      Facet_pt[3]->set_vertex_pt(2, Vertex_pt[5]);

      // icosa_facet[3].resize(3);
      // icosa_facet[3][0] = 0;
      // icosa_facet[3][1] = 4;
      // icosa_facet[3][2] = 5;

      Facet_pt[4] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[4]->set_vertex_pt(0, Vertex_pt[0]);
      Facet_pt[4]->set_vertex_pt(1, Vertex_pt[8]);
      Facet_pt[4]->set_vertex_pt(2, Vertex_pt[4]);

      // icosa_facet[4].resize(3);
      // icosa_facet[4][0] = 0;
      // icosa_facet[4][1] = 8;
      // icosa_facet[4][2] = 4;

      Facet_pt[5] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[5]->set_vertex_pt(0, Vertex_pt[5]);
      Facet_pt[5]->set_vertex_pt(1, Vertex_pt[11]);
      Facet_pt[5]->set_vertex_pt(2, Vertex_pt[10]);

      // icosa_facet[5].resize(3);
      // icosa_facet[5][0] = 5;
      // icosa_facet[5][1] = 11;
      // icosa_facet[5][2] = 10;

      Facet_pt[6] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[6]->set_vertex_pt(0, Vertex_pt[5]);
      Facet_pt[6]->set_vertex_pt(1, Vertex_pt[2]);
      Facet_pt[6]->set_vertex_pt(2, Vertex_pt[11]);

      // icosa_facet[6].resize(3);
      // icosa_facet[6][0] = 5;
      // icosa_facet[6][1] = 2;
      // icosa_facet[6][2] = 11;

      Facet_pt[7] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[7]->set_vertex_pt(0, Vertex_pt[4]);
      Facet_pt[7]->set_vertex_pt(1, Vertex_pt[2]);
      Facet_pt[7]->set_vertex_pt(2, Vertex_pt[5]);

      // icosa_facet[7].resize(3);
      // icosa_facet[7][0] = 4;
      // icosa_facet[7][1] = 2;
      // icosa_facet[7][2] = 5;

      Facet_pt[8] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[8]->set_vertex_pt(0, Vertex_pt[4]);
      Facet_pt[8]->set_vertex_pt(1, Vertex_pt[9]);
      Facet_pt[8]->set_vertex_pt(2, Vertex_pt[2]);

      // icosa_facet[8].resize(3);
      // icosa_facet[8][0] = 4;
      // icosa_facet[8][1] = 9;
      // icosa_facet[8][2] = 2;

      Facet_pt[9] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[9]->set_vertex_pt(0, Vertex_pt[8]);
      Facet_pt[9]->set_vertex_pt(1, Vertex_pt[9]);
      Facet_pt[9]->set_vertex_pt(2, Vertex_pt[4]);

      // icosa_facet[9].resize(3);
      // icosa_facet[9][0] = 8;
      // icosa_facet[9][1] = 9;
      // icosa_facet[9][2] = 4;

      Facet_pt[10] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[10]->set_vertex_pt(0, Vertex_pt[6]);
      Facet_pt[10]->set_vertex_pt(1, Vertex_pt[9]);
      Facet_pt[10]->set_vertex_pt(2, Vertex_pt[8]);

      // icosa_facet[10].resize(3);
      // icosa_facet[10][0] = 6;
      // icosa_facet[10][1] = 9;
      // icosa_facet[10][2] = 8;

      Facet_pt[11] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[11]->set_vertex_pt(0, Vertex_pt[1]);
      Facet_pt[11]->set_vertex_pt(1, Vertex_pt[6]);
      Facet_pt[11]->set_vertex_pt(2, Vertex_pt[8]);

      // icosa_facet[11].resize(3);
      // icosa_facet[11][0] = 1;
      // icosa_facet[11][1] = 6;
      // icosa_facet[11][2] = 8;

      Facet_pt[12] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[12]->set_vertex_pt(0, Vertex_pt[1]);
      Facet_pt[12]->set_vertex_pt(1, Vertex_pt[7]);
      Facet_pt[12]->set_vertex_pt(2, Vertex_pt[6]);

      // icosa_facet[12].resize(3);
      // icosa_facet[12][0] = 1;
      // icosa_facet[12][1] = 7;
      // icosa_facet[12][2] = 6;

      Facet_pt[13] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[13]->set_vertex_pt(0, Vertex_pt[10]);
      Facet_pt[13]->set_vertex_pt(1, Vertex_pt[7]);
      Facet_pt[13]->set_vertex_pt(2, Vertex_pt[1]);

      // icosa_facet[13].resize(3);
      // icosa_facet[13][0] = 10;
      // icosa_facet[13][1] = 7;
      // icosa_facet[13][2] = 1;

      Facet_pt[14] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[14]->set_vertex_pt(0, Vertex_pt[10]);
      Facet_pt[14]->set_vertex_pt(1, Vertex_pt[11]);
      Facet_pt[14]->set_vertex_pt(2, Vertex_pt[7]);

      // icosa_facet[14].resize(3);
      // icosa_facet[14][0] = 10;
      // icosa_facet[14][1] = 11;
      // icosa_facet[14][2] = 7;

      Facet_pt[15] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[15]->set_vertex_pt(0, Vertex_pt[11]);
      Facet_pt[15]->set_vertex_pt(1, Vertex_pt[3]);
      Facet_pt[15]->set_vertex_pt(2, Vertex_pt[7]);

      // icosa_facet[15].resize(3);
      // icosa_facet[15][0] = 11;
      // icosa_facet[15][1] = 3;
      // icosa_facet[15][2] = 7;

      Facet_pt[16] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[16]->set_vertex_pt(0, Vertex_pt[7]);
      Facet_pt[16]->set_vertex_pt(1, Vertex_pt[3]);
      Facet_pt[16]->set_vertex_pt(2, Vertex_pt[6]);

      // icosa_facet[16].resize(3);
      // icosa_facet[16][0] = 7;
      // icosa_facet[16][1] = 3;
      // icosa_facet[16][2] = 6;

      Facet_pt[17] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[17]->set_vertex_pt(0, Vertex_pt[6]);
      Facet_pt[17]->set_vertex_pt(1, Vertex_pt[3]);
      Facet_pt[17]->set_vertex_pt(2, Vertex_pt[9]);

      // icosa_facet[17].resize(3);
      // icosa_facet[17][0] = 6;
      // icosa_facet[17][1] = 3;
      // icosa_facet[17][2] = 9;

      Facet_pt[18] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[18]->set_vertex_pt(0, Vertex_pt[9]);
      Facet_pt[18]->set_vertex_pt(1, Vertex_pt[3]);
      Facet_pt[18]->set_vertex_pt(2, Vertex_pt[2]);

      // icosa_facet[18].resize(3);
      // icosa_facet[18][0] = 9;
      // icosa_facet[18][1] = 3;
      // icosa_facet[18][2] = 2;

      Facet_pt[19] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[19]->set_vertex_pt(0, Vertex_pt[2]);
      Facet_pt[19]->set_vertex_pt(1, Vertex_pt[3]);
      Facet_pt[19]->set_vertex_pt(2, Vertex_pt[11]);

      // icosa_facet[19].resize(3);
      // icosa_facet[19][0] = 2;
      // icosa_facet[19][1] = 3;
      // icosa_facet[19][2] = 11;

      // Set one-based boundary IDs
      unsigned one_based_boundary_id = 1;
      for (unsigned f = 0; f < n_facet; f++)
      {
        Facet_pt[f]->set_one_based_boundary_id(one_based_boundary_id);
      }

      // Identify point in hole
      Vector<double> inner_point(3, 0.0);
      set_hole_for_tetgen(inner_point);
    }
  };

} // namespace oomph
