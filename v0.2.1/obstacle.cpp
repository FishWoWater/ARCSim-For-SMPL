/*
  Copyright ©2013 The Regents of the University of California
  (Regents). All Rights Reserved. Permission to use, copy, modify, and
  distribute this software and its documentation for educational,
  research, and not-for-profit purposes, without fee and without a
  signed licensing agreement, is hereby granted, provided that the
  above copyright notice, this paragraph and the following two
  paragraphs appear in all copies, modifications, and
  distributions. Contact The Office of Technology Licensing, UC
  Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620,
  (510) 643-7201, for commercial licensing opportunities.

  IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT,
  INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
  LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS
  DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF THE POSSIBILITY
  OF SUCH DAMAGE.

  REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING
  DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS
  IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
*/

#include "io.hpp"
#include "obstacle.hpp"
#include "util.hpp"
#include <cstdio>

using namespace std;

Mesh& Obstacle::get_mesh() {
    return curr_state_mesh;
}

const Mesh& Obstacle::get_mesh() const {
    return curr_state_mesh;
}

Mesh& Obstacle::get_mesh(double time) {
    if (time > end_time)
        delete_mesh(curr_state_mesh);
    if (time < start_time || time > end_time)
        return curr_state_mesh;
    if (!activated)
        curr_state_mesh = deep_copy(base_mesh);
    if (transform_spline) {
        DTransformation dtrans = get_dtrans(*transform_spline, time);
        Mesh &mesh = curr_state_mesh;
        for (int n = 0; n < curr_state_mesh.nodes.size(); n++)
            mesh.nodes[n]->x = apply_dtrans(dtrans, base_mesh.nodes[n]->x,
                                            &mesh.nodes[n]->v);
        compute_ws_data(mesh);
    }
    if (!activated)
        update_x0(curr_state_mesh);
    activated = true;
    return curr_state_mesh;
}


Mesh& Obstacle::get_smpl_mesh(double time) {
    // 1. interpolate the quaternion in transformation 
    // 2. convert quaternion to axis-angle 
    // 3. parse axis-angle poses and betas to smpl model;
    // std::cout << "entering get smpl mesh" << std::endl;
    Quaternion x; 
    // x.to_axisangle; 
    if (time > end_time)
        delete_mesh(curr_state_mesh);
    if (time < start_time || time > end_time)
        return curr_state_mesh;
    if (!activated)
        curr_state_mesh = deep_copy(base_mesh);
    if (smpl_motion) {
//        std::vector<Quaternion> smpl_quaternions = get_smpl_quaternions(*smpl_motion, time);
        Transformation smpl_transformation = get_smpl_transformation(*smpl_motion, time);
        std::vector<Quaternion> smpl_quaternions = smpl_transformation.rotations;
        Vec3 smpl_translation = smpl_transformation.translation;
        // used for interpolation
        std::vector<double> dynamic_betas = smpl_transformation.dynamic_betas;
        // std::cout << "smpl translation:" << smpl_translation << std::endl;
        Eigen::VectorXf thetas, dyna_betas;
        thetas.resize(72);
        dyna_betas.resize(10);
        for(int i=0;i<10;i++){
            dyna_betas(i) = dynamic_betas[i];
        }
        // parse thetas from the quaternions
        for(int i=0;i<smpl_quaternions.size();i++){
            auto aa_pair = smpl_quaternions[i].to_axisangle();
            Vec3 axis = aa_pair.first; 
            double angle = aa_pair.second;
            thetas[3 * i] = axis[0] * angle;
            thetas[3 * i + 1] = axis[1] * angle; 
            thetas[3 * i + 2] = axis[2] * angle;
        }
        // forward kinematics 
        smpl->setAllPoses(thetas);

        smpl->setAllShapes(dyna_betas);
        smpl->updateModel();
//         std::cout << "smpl model updated" << std::endl;
//        smpl->saveToOBJ("rest_smpl.obj");
//        exit(0);

        Eigen::MatrixXf vertices = smpl->mVTemp2;
        // assign the vertices to mesh 
        Mesh &mesh = curr_state_mesh;
        for (int n = 0; n < curr_state_mesh.nodes.size(); n++) {
            // add transformation here
//            mesh.nodes[n]->x[0] = vertices(n, 0);
//            mesh.nodes[n]->x[1] = vertices(n, 1);
//            mesh.nodes[n]->x[2] = vertices(n, 2);
            mesh.nodes[n]->x[0] = vertices(n, 0) + smpl_translation[0];
            mesh.nodes[n]->x[1] = vertices(n, 1) + smpl_translation[1];
            mesh.nodes[n]->x[2] = vertices(n, 2) + smpl_translation[2];
        }
        // std::cout << "smpl mesh assigned to arcsim mesh" << std::endl; 
        compute_ws_data(mesh);
    }
    if (!activated)
        update_x0(curr_state_mesh);
    activated = true;
    // std::cout << "smpl mesh updated" << std::endl; 
    return curr_state_mesh;
}

void Obstacle::blend_with_previous (double t, double dt, double blend) {
    const Motion *spline = transform_spline;
    Transformation trans = (spline)
                         ? get_trans(*spline, t)
                           * inverse(get_trans(*spline, t-dt))
                         : identity_();
    Mesh &mesh = curr_state_mesh;
    for (int n = 0; n < mesh.nodes.size(); n++) {
        Node *node = mesh.nodes[n];
        Vec3 x0 = trans.apply(node->x0);
        node->x = x0 + blend*(node->x - x0);
    }
    compute_ws_data(mesh);
}

void Obstacle::smpl_blend_with_previous (double t, double dt, double blend) {
    const Motion *spline = smpl_motion;
//    Transformation trans = (spline)
//                           ? get_trans(*spline, t)
//                             * inverse(get_trans(*spline, t-dt))
//                           : identity_();
    Mesh &mesh = curr_state_mesh;
    for (int n = 0; n < mesh.nodes.size(); n++) {
        Node *node = mesh.nodes[n];
//        Vec3 x0 = trans.apply(node->x0);
        Vec3 x0 = node -> x0;
        node->x = x0 + blend*(node->x - x0);
    }
    compute_ws_data(mesh);
}

