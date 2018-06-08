// Iterative solver for kinematic moves
//
// Copyright (C) 2018  Kevin O'Connor <kevin@koconnor.net>
//
// This file may be distributed under the terms of the GNU GPLv3 license.

#include <math.h> // sqrt
#include <stdlib.h> // malloc
#include <string.h> // memset
#include "compiler.h" // __visible
#include "itersolve.h" // struct coord
#include "pyhelper.h" // errorf
#include "stepcompress.h" // queue_append_start


/****************************************************************
 * Kinematic moves
 ****************************************************************/

struct move {
    double print_time, move_t;
    double accel_t, cruise_t;
    double start_v, cruise_v;
    double half_accel;
    double cruise_start_d, decel_start_d;
    struct coord start_pos, axes_r;
};

struct move * __visible
move_alloc(void)
{
    struct move *m = malloc(sizeof(*m));
    memset(m, 0, sizeof(*m));
    return m;
}

// Populate a 'struct move' with a velocity trapezoid
void __visible
move_fill(struct move *m, double print_time
          , double accel_t, double cruise_t, double decel_t
          , double start_pos_x, double start_pos_y, double start_pos_z
          , double axes_d_x, double axes_d_y, double axes_d_z
          , double start_v, double cruise_v, double accel)
{
    m->print_time = print_time;
    m->move_t = accel_t + cruise_t + decel_t;
    m->accel_t = accel_t;
    m->cruise_t = cruise_t;
    m->start_v = start_v;
    m->cruise_v = cruise_v;
    m->half_accel = .5 * accel;
    m->cruise_start_d = accel_t * .5 * (cruise_v + start_v);
    m->decel_start_d = m->cruise_start_d + cruise_t * cruise_v;
    m->start_pos.x = start_pos_x;
    m->start_pos.y = start_pos_y;
    m->start_pos.z = start_pos_z;
    double inv_move_d = 1. / sqrt(axes_d_x*axes_d_x + axes_d_y*axes_d_y
                                  + axes_d_z*axes_d_z);
    m->axes_r.x = axes_d_x * inv_move_d;
    m->axes_r.y = axes_d_y * inv_move_d;
    m->axes_r.z = axes_d_z * inv_move_d;
}

// Return the distance moved given a time in a move
static double
move_get_distance(struct move *m, double move_time)
{
    if (move_time < m->accel_t)
        return (m->start_v + move_time*m->half_accel) * move_time;
    move_time -= m->accel_t;
    if (move_time < m->cruise_t)
        return m->cruise_start_d + m->cruise_v * move_time;
    move_time -= m->cruise_t;
    double avg_v = m->cruise_v - move_time*m->half_accel;
    return m->decel_start_d + avg_v * move_time;
}

// Return the XYZ coordinates given a time in a move
inline struct coord
move_get_coord(struct move *m, double move_time)
{
    double move_dist = move_get_distance(m, move_time);
    return (struct coord) {
        .x = m->start_pos.x + m->axes_r.x * move_dist,
        .y = m->start_pos.y + m->axes_r.y * move_dist,
        .z = m->start_pos.z + m->axes_r.z * move_dist };
}


/****************************************************************
 * Iterative solver
 ****************************************************************/

int32_t __visible
itersolve_gen_steps(struct stepper_kinematics *sk, struct move *m)
{
    struct stepcompress *sc = sk->sc;
    sk_callback calc_position = sk->calc_position;
    double half_step_dist = .5 * sk->step_dist;
    double mcu_freq = stepcompress_get_mcu_freq(sc);
    double last_pos = sk->commanded_pos;
    double low_guess_pos = last_pos, high_guess_pos = last_pos;
    double last_time = 0., low_guess_time = 0., high_guess_time = 0.;
    double seek_time_delta = 0.000100;
    int sdir = stepcompress_get_step_dir(sc);
    int steps = 0;
    struct queue_append qa = queue_append_start(sc, m->print_time, .5);
    for (;;) {
        // Determine if next step is in forward or reverse direction
        double dist = high_guess_pos - last_pos;
        if (high_guess_time <= low_guess_time || fabs(dist) < half_step_dist) {
            if (high_guess_time >= m->move_t)
                // At end of move
                break;
            // Need to increase high_guess_time
            low_guess_time = high_guess_time;
            low_guess_pos = high_guess_pos;
            high_guess_time = last_time + seek_time_delta;
            seek_time_delta += seek_time_delta;
            if (high_guess_time > m->move_t)
                high_guess_time = m->move_t;
            high_guess_pos = calc_position(sk, m, high_guess_time);
            continue;
        }
        int next_sdir = dist > 0.;
        if (next_sdir != sdir) {
            int ret = queue_append_set_next_step_dir(&qa, next_sdir);
            if (ret)
                return ret;
            sdir = next_sdir;
        }
        // Find step using "false position" method
        double target = last_pos + (sdir ? half_step_dist : -half_step_dist);
        double min_guess_time = low_guess_time;
        double min_guess_dpos = low_guess_pos - target;
        double max_guess_time = high_guess_time, guess_time = high_guess_time;
        double max_guess_dpos = high_guess_pos - target;
        double guess_pos = high_guess_pos;
        int max_guess_sign = signbit(max_guess_dpos);
        for (;;) {
            double minmax_delta = max_guess_dpos - min_guess_dpos;
            if (!minmax_delta)
                break;
            double next_guess_time = ((min_guess_time*max_guess_dpos
                                       - max_guess_time*min_guess_dpos)
                                      / minmax_delta);
            if (fabs(next_guess_time - guess_time) <= .000000001)
                break;
            guess_time = next_guess_time;
            guess_pos = calc_position(sk, m, guess_time);
            double guess_dpos = guess_pos - target;
            if (!guess_dpos)
                break;
            int guess_sign = signbit(guess_dpos);
            if (guess_sign == max_guess_sign) {
                max_guess_time = guess_time;
                max_guess_dpos = guess_dpos;
            } else {
                min_guess_time = guess_time;
                min_guess_dpos = guess_dpos;
            }
        }
        // Add step at given time
        int ret = queue_append(&qa, guess_time * mcu_freq);
        if (ret)
            return ret;
        last_pos = sdir ? target + half_step_dist : target - half_step_dist;
        steps += sdir ? 1 : -1;
        low_guess_pos = guess_pos;
        seek_time_delta = guess_time - last_time;
        if (seek_time_delta < .000000001)
            seek_time_delta = .000000001;
        last_time = low_guess_time = guess_time;
    }
    queue_append_finish(qa);
    sk->commanded_pos = last_pos;
    return steps;
}

void __visible
itersolve_set_stepcompress(struct stepper_kinematics *sk
                           , struct stepcompress *sc, double step_dist)
{
    sk->sc = sc;
    sk->step_dist = step_dist;
}

void __visible
itersolve_set_commanded_pos(struct stepper_kinematics *sk, double pos)
{
    sk->commanded_pos = pos;
}
