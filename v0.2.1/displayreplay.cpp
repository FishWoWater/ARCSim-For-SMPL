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

#include "displayreplay.hpp"

#include "hylc/hylc_conf.hpp"
#include "hylc/hylc.hpp"

#include "conf.hpp"
#include "display.hpp"
#include "io.hpp"
#include "misc.hpp"
#include "opengl.hpp"
#include <cstdio>
#include <fstream>
#include <boost/filesystem.hpp>
using namespace std;

#ifndef NO_OPENGL

static string inprefix, outprefix;
static int frameskip;

static bool running = true;

static void reload () {
    int fullframe = ::frame*::frameskip;
    sim.time = fullframe * sim.frame_time;
    // std::cout << "cloth file" << stringf("%s/%04d_00.obj",inprefix.c_str(), fullframe) << std::endl;
    if(!boost::filesystem::exists(stringf("%s/%04d_00.obj",inprefix.c_str(), fullframe))){
        std::cout << "can not find cloth file" << std::endl;
        if (::frame == 0)
            exit(EXIT_FAILURE);
        if (!outprefix.empty())
            exit(EXIT_SUCCESS);
        ::frame = 0;
        reload();
    }
    load_objs(sim.cloth_meshes, stringf("%s/%04d",inprefix.c_str(), fullframe));
    if (sim.cloth_meshes[0]->verts.empty() and !sim.is_smpl) {
        if (::frame == 0)
            exit(EXIT_FAILURE);
        if (!outprefix.empty())
            exit(EXIT_SUCCESS);
        ::frame = 0;
        reload();
    }
    for (int o = 0; o < sim.obstacles.size(); o++) {
        if(!sim.is_smpl)    sim.obstacles[o].get_mesh(sim.time);
        else{
//            double decay_time = 0.1, blend = sim.step_time / decay_time;
//            blend = blend / (1 + blend);
            sim.obstacles[o].get_smpl_mesh(sim.time);
//            sim.obstacles[o].smpl_blend_with_previous(sim.time, sim.step_time, blend);
        }
    }
}

static void idle () {
    if (!running)
        return;
    fps.tick();
    if (!outprefix.empty()) {
        char filename[256];
        snprintf(filename, 256, "%s/%04d.png", outprefix.c_str(), ::frame);
        save_screenshot(filename);
    }
    ::frame++;
    reload();
    fps.tock();
    redisplay();
}

static void keyboard (unsigned char key, int x, int y) {
    unsigned char esc = 27, space = ' ';
    if (key == esc) {
        exit(0);
    } else if (key == space) {
        running = !running;
    } else if (key == 'p') {
        std::string filename = "TMP_c";
        int c = 0;
        for (auto & cloth : sim.cloths) {
            hylc::hylc_write_strains(filename + std::to_string(c) + ".txt", cloth);
            c++;
        }
    }
}

static void special (int key, int x, int y) {
    bool shift = glutGetModifiers() & GLUT_ACTIVE_SHIFT,
         alt = glutGetModifiers() & GLUT_ACTIVE_ALT;
    int delta = alt ? 100 : shift ? 10 : 1;
    if (key == GLUT_KEY_LEFT) {
        ::frame -= delta;
        reload();
    } else if (key == GLUT_KEY_RIGHT) {
        ::frame += delta;
        reload();
    } else if (key == GLUT_KEY_HOME) {
        ::frame = 0;
        reload();
    }
    redisplay();
}

void display_replay (const vector<string> &args) {
    if (args.size() < 1 || args.size() > 2) {
        cout << "Replays the results of a simulation." << endl;
        cout << "Arguments:" << endl;
        cout << "    <out-dir>: Directory containing simulation output files"
             << endl;
        cout << "    <sshot-dir> (optional): Directory to save images" << endl;
        exit(EXIT_FAILURE);
    }
    ::inprefix = args[0];
    ::outprefix = args.size()>1 ? args[1] : "";
    ::frameskip = 1;
    if (!::outprefix.empty())
        ensure_existing_directory(::outprefix);
    char config_backup_name[256];
    snprintf(config_backup_name, 256, "%s/%s", inprefix.c_str(), "conf.json");
    load_json(config_backup_name, sim);
    prepare(sim);
    reload();
    GlutCallbacks cb;
    cb.idle = idle;
    cb.keyboard = keyboard;
    cb.special = special;
    run_glut(cb);
}


#else

void display_replay (const vector<string> &args) {opengl_fail();}

#endif // NO_OPENGL
