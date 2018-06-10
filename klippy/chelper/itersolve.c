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

struct timepos {
    double time, position;
};

// Find step using "false position" method
static struct timepos
itersolve_find_step(struct stepper_kinematics *sk, struct move *m
                    , struct timepos low, struct timepos high
                    , double target)
{
    sk_callback calc_position = sk->calc_position;
    struct timepos best_guess = high;
    low.position -= target;
    high.position -= target;
    int high_sign = signbit(high.position);
    for (;;) {
        double range_delta = high.position - low.position;
        if (!range_delta)
            break;
        double guess_time = ((low.time*high.position - high.time*low.position)
                             / range_delta);
        if (fabs(guess_time - best_guess.time) <= .000000001)
            break;
        best_guess.time = guess_time;
        best_guess.position = calc_position(sk, m, guess_time);
        double guess_position = best_guess.position - target;
        if (!guess_position)
            break;
        int guess_sign = signbit(guess_position);
        if (guess_sign == high_sign) {
            high.time = guess_time;
            high.position = guess_position;
        } else {
            low.time = guess_time;
            low.position = guess_position;
        }
    }
    return best_guess;
}

int32_t __visible
itersolve_gen_steps(struct stepper_kinematics *sk, struct move *m)
{
    struct stepcompress *sc = sk->sc;
    sk_callback calc_position = sk->calc_position;
    double half_step = .5 * sk->step_dist;
    double mcu_freq = stepcompress_get_mcu_freq(sc);
    struct timepos last = { 0., sk->commanded_pos }, low = last, high = last;
    double seek_time_delta = 0.000100;
    int sdir = stepcompress_get_step_dir(sc);
    int steps = 0;
    struct queue_append qa = queue_append_start(sc, m->print_time, .5);
    for (;;) {
        // Determine if next step is in forward or reverse direction
        double dist = high.position - last.position;
        if (high.time <= low.time || fabs(dist) < half_step) {
            if (high.time >= m->move_t)
                // At end of move
                break;
            // Need to increase next step search range
            low = high;
            high.time = last.time + seek_time_delta;
            seek_time_delta += seek_time_delta;
            if (high.time > m->move_t)
                high.time = m->move_t;
            high.position = calc_position(sk, m, high.time);
            continue;
        }
        int next_sdir = dist > 0.;
        if (next_sdir != sdir) {
            int ret = queue_append_set_next_step_dir(&qa, next_sdir);
            if (ret)
                return ret;
            sdir = next_sdir;
        }
        // Find step
        double target = last.position + (sdir ? half_step : -half_step);
        struct timepos next = itersolve_find_step(sk, m, low, high, target);
        // Add step at given time
        int ret = queue_append(&qa, next.time * mcu_freq);
        if (ret)
            return ret;
        steps += sdir ? 1 : -1;
        seek_time_delta = next.time - last.time;
        if (seek_time_delta < .000000001)
            seek_time_delta = .000000001;
        last.position = target + (sdir ? half_step : -half_step);
        last.time = next.time;
        low = next;
    }
    queue_append_finish(qa);
    sk->commanded_pos = last.position;
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
