# TMC2130 configuration
#
# Copyright (C) 2018  Kevin O'Connor <kevin@koconnor.net>
#
# This file may be distributed under the terms of the GNU GPLv3 license.
import math
import logging

TMC_FREQUENCY=13200000.
GCONF_EN_PWM_MODE=1<<2
GCONF_DIAG1_STALL=1<<8
REG_GCONF=0x00
REG_TCOOLTHRS=0x14
REG_COOLCONF=0x6d
REG_PWMCONF=0x70
REG_CHOPCONF=0x6c
REG_IHOLD_IRUN=0x10
REG_TPOWERDOWN=0x11
REG_TPWMTHRS=0x13
REG_MSLUTSTART=0x69
REG_MSLUT0=0x60
REG_MSLUTSEL = 0x68

# Constants for sine wave correction
TMC_WAVE_FACTOR_MIN = 1.005
TMC_WAVE_FACTOR_MAX = 1.3
TMC_WAVE_AMP = 247

class TMC2130:
    def __init__(self, config):
        self.printer = config.get_printer()
        self.stepper_name = config.get_name().split()[1]
        # Pin setup
        ppins = self.printer.lookup_object("pins")
        cs_pin = config.get('cs_pin')
        cs_pin_params = ppins.lookup_pin('digital_out', cs_pin)
        if cs_pin_params['invert']:
            raise pins.error("tmc2130 can not invert pin")
        self.mcu = cs_pin_params['chip']
        pin = cs_pin_params['pin']
        self.oid = self.mcu.create_oid()
        self.mcu.add_config_cmd(
            "config_spi oid=%d bus=%d pin=%s mode=%d rate=%d shutdown_msg=" % (
                self.oid, 0, cs_pin_params['pin'], 3, 4000000))
		# Calculate current	and populate IHOLD_IRUN			
        run_current = config.getfloat('run_current', above=0., maxval=2.)
        hold_current = config.getfloat('hold_current', run_current,
                                       above=0., maxval=2.)'
        sense_resistor = config.getfloat('sense_resistor', 0.110, above=0.)									   
        vsense = False
        irun = self.current_bits(run_current, sense_resistor, vsense)
        ihold = self.current_bits(hold_current, sense_resistor, vsense)
        if irun < 16 and ihold < 16:
            vsense = True
            irun = self.current_bits(run_current, sense_resistor, vsense)
            ihold = self.current_bits(hold_current, sense_resistor, vsense)
        iholddelay = config.getint('driver_IHOLDDELAY', 8, minval=0, maxval=15)
        tpowerdown = config.getint('driver_TPOWERDOWN', 0, minval=0, maxval=255)
        # Populate GCONF
        en_pwm_mode = config.getboolean('stealthchop', False)
		# Populate TPWMTHRS	
        sc_velocity = config.getfloat('stealthchop_threshold', 0., minval=0.)
        sc_threshold = self.velocity_to_clock(config, sc_velocity)
        # Populate TCOOLTHRS
        coolstep_velocity = config.getfloat('coolstep_threshold', 0., minval=0.)
        coolstep_threshold = self.velocity_to_clock(config, coolstep_velocity) 
        # Populate CHOPCONF
        interpolate = config.getboolean('interpolate', True)
        steps = {'256': 0, '128': 1, '64': 2, '32': 3, '16': 4,
                 '8': 5, '4': 6, '2': 7, '1': 8}
        self.mres = config.getchoice('microsteps', steps)
        pwm_sync = config.getint('driver_PWM_SYNC', 0, minval=1, maxval=15)       
        blank_time_select = config.getint('driver_BLANK_TIME_SELECT', 1,
                                          minval=0, maxval=3)
        chm = config.getint('driver_CHOPPER_MODE', 0, minval=0, maxval=1)
        disfdcc = config.getboolean('driver_FAST_DECAY', False)
        rndtf = config.getboolean('driver_RANDOM_OFF_TIME', False)
        hend = config.getint('driver_HEND', 7, minval=0, maxval=15)
        hstrt = config.getint('driver_HSTRT', 0, minval=0, maxval=7)
        toff = config.getint('driver_TOFF', 4, minval=1, maxval=15)
        # Populate COOLCONF
        sfilt = config.getint('driver_SG_FILTER', 0, minval=0, maxval=1)
        sgt = config.getint('driver_SGT', 0, minval=-64, maxval=63) & 0x7f
        seimin = config.getint('driver_CURRENT_MIN', 0, minval=0, maxval= 1)
        sedn = config.getint('driver_CURRENT_SPEED', 0, minval=0, maxval=3)
        semax = config.getint('driver_SGT_HYST', 0, minval=0, maxval=15)
        seup = config.getint('driver_CURRENT_STEP', 1, minval=0, maxval=3)
        semin = config.getint('driver_SGT_MIN', 0, minval=0, maxval=15)
        # Populate PWMCONF
        freewheel = config.getint('driver_STANDSTILL_MODE', 0, minval=0, maxval=3)
        pwm_sym = config.getboolean('driver_SYMMETRIC_PWM', False)
        pwm_scale = config.getboolean('driver_PWM_AUTOSCALE', True)
        pwm_freq = config.getint('driver_PWM_FREQ', 1, minval=0, maxval=3)
        pwm_grad = config.getint('driver_PWM_GRAD', 4, minval=0, maxval=255)
        pwm_ampl = config.getint('driver_PWM_AMPL', 128, minval=0, maxval=255)
        # configure GCONF
        self.add_config_cmd(REG_GCONF, en_pwm_mode << 2)
        # configure IHOLD_IRUN
        self.add_config_cmd(REG_IHOLD_IRUN, ihold | (irun << 8) | (iholddelay << 16))
        # configure TPOWERDOWN
        self.add_config_cmd(REG_TPOWERDOWN, tpowerdown)
        # configure TPWMTHRS
        self.add_config_cmd(REG_TPWMTHRS, max(0, min(0xfffff, sc_threshold)))
		# configure TCOOLTHRS
        self.add_config_cmd(REG_TCOOLTHRS, coolstep_threshold)
        # configure CHOPCONF
        self.add_config_cmd(
            REG_CHOPCONF, toff | (hstrt << 4) | (hend << 7) | (disfdcc << 12) | (rndtf << 13) 
            | (chm << 14) | (blank_time_select << 15) | (vsense << 17) | (pwm_sync << 20) 
            | (self.mres << 24) | (interpolate << 28))
        # configure COOLCONF
        self.add_config_cmd(
            REG_COOLCONF, semin | (seup << 5) | (semax << 8) | (sedn << 13) | (seimin << 15) 
            | (sgt << 16) | (sfilt << 24))
        # configure PWMCONF
        self.add_config_cmd(
            REG_PWMCONF, pwm_ampl | (pwm_grad << 8) | (pwm_freq << 16) | (pwm_scale << 18) 
            |(pwm_sym << 19) | (freewheel << 20))
        # configure Linearity Correction
        wave_factor = config.getfloat('linearity_correction', 0., minval=TMC_WAVE_FACTOR_MIN, maxval=TMC_WAVE_FACTOR_MAX)
        self.set_wave(wave_factor, True)
        # Linearity Correction GCODE setup
        self.gcode = self.printer.lookup_object('gcode')
        self.gcode.register_mux_command("TMC_SET_WAVE", "STEPPER", self.stepper_name,
                                        self.cmd_TMC_SET_WAVE,
                                        desc=self.cmd_TMC_SET_WAVE_help)
        # Allow virtual endstop to be created
        self.diag1_pin = config.get('diag1_pin', None)
        ppins.register_chip("_".join(config.get_name().split()[:2]), self)
        self.send_spi_cmd = None
        self.mcu.add_config_object(self)
    def add_config_cmd(self, addr, val):
        self.mcu.add_config_cmd("spi_send oid=%d data=%02x%08x" % (
            self.oid, (addr | 0x80) & 0xff, val & 0xffffffff), is_init=True)
    def current_bits(self, current, sense_resistor, vsense_on):
        sense_resistor += 0.020
        vsense = 0.32
        if vsense_on:
            vsense = 0.18
        cs = int(32. * current * sense_resistor * math.sqrt(2.) / vsense
                 - 1. + .5)
        return max(0, min(31, cs))
    def velocity_to_clock(self, config, velocity):
        if not velocity:
            return 0
        stepper_name = config.get_name().split()[1]
        stepper_config = config.getsection(stepper_name)
        step_dist = stepper_config.getfloat('step_distance')
        step_dist_256 = step_dist / (1 << self.mres)
        return int(TMC_FREQUENCY * step_dist_256 / velocity + .5)
    def setup_pin(self, pin_params):
        if (pin_params['pin'] != 'virtual_endstop'
            or pin_params['type'] != 'endstop'):
            raise pins.error("tmc2130 virtual endstop only useful as endstop")
        if pin_params['invert'] or pin_params['pullup']:
            raise pins.error("Can not pullup/invert tmc2130 virtual endstop")
        return TMC2130VirtualEndstop(self)
    def build_config(self):
        cmd_queue = self.mcu.alloc_command_queue()
        self.spi_send_cmd = self.mcu.lookup_command(
            "spi_send oid=%c data=%*s", cq=cmd_queue)
    def set_register(self, addr, val):
        data = [(addr | 0x80) & 0xff, (val >> 24) & 0xff, (val >> 16) & 0xff,
                (val >> 8) & 0xff, val & 0xff]
        self.spi_send_cmd.send([self.oid, data])
    def set_wave(self, fac, init=False):
        if fac < TMC_WAVE_FACTOR_MIN:
             fac = 0.0
        elif fac > TMC_WAVE_FACTOR_MAX:
            fac = TMC_WAVE_FACTOR_MAX
        error = None
        vA = 0
        prevA = 0
        delta0 = 0
        delta1 = 1
        w = [1, 1, 1, 1]
        x = [255, 255, 255]
        seg = 0
        bitVal = 0
        deltaA = 0
        reg = 0
        # configure MSLUTSTART
        if init:
            self.add_config_cmd(REG_MSLUTSTART, (TMC_WAVE_AMP << 16))
        else:
            self.set_register(REG_MSLUTSTART, (TMC_WAVE_AMP << 16))
        for i in range(256):
            if (i & 31) == 0:
                reg = 0
            if fac == 0.:
                # default TMC wave
                vA = int((TMC_WAVE_AMP + 1) * math.sin((2*math.pi*i + math.pi)/1024) + .5) - 1
            else:
                # corrected wave
                vA = int(TMC_WAVE_AMP * math.pow(math.sin(2*math.pi*i/1024), fac) + .5)
            deltaA = vA - prevA
            prevA = vA
            bitVal = -1
            if deltaA == delta0:
                bitVal = 0
            elif deltaA == delta1:
                bitVal = 1
            else:
                if deltaA < delta0:
                    #switch w bit down
                    bitVal = 0
                    if deltaA == -1:
                        delta0 = -1
                        delta1 = 0
                        w[seg+1] = 0
                    elif deltaA == 0:
                        delta0 = 0
                        delta1 = 1
                        w[seg+1] = 1
                    elif deltaA == 1:
                        delta0 = 1
                        delta1 = 2
                        w[seg+1] = 2
                    else:
                        bitVal = -1
                    if bitVal >= 0:
                        x[seg] = i
                        seg += 1
                elif deltaA > delta1:
                    #switch w bit up
                    bitVal = 1
                    if deltaA == 1:
                        delta0 = 0
                        delta1 = 1
                        w[seg+1] = 1
                    elif deltaA == 2:
                        delta0 = 1
                        delta1 = 2
                        w[seg+1] = 2
                    elif deltaA == 3:
                        delta0 = 2
                        delta1 = 3
                        w[seg+1] = 3
                    else:
                        bitVal = -1
                    if bitVal >= 0:
                        x[seg] = i
                        seg += 1
            if bitVal < 0:
                # delta out of range
                error = "TMC2130: Error setting Sine Wave, Delta Out of Range"
                break
            if seg > 3: # TODO: should this be greater than 2?
                # segment out of range
                error = "TMC2130: Error setting Sine Wave, Segment Out of Range"
                break
            if bitVal == 1:
                reg |= 0x80000000
            if (i & 31) == 31:
                #configure MSLUT
                if init:
                    self.add_config_cmd(REG_MSLUT0 + ((i >> 5) & 7), reg)
                else:
                    self.set_register(REG_MSLUT0 + ((i >> 5) & 7), reg)
            else:
                reg >>= 1
        success_msg = "TMC2130: Wave factor on stepper [%s] set to: %f" % (self.stepper_name, fac)
        # configure MSLUTSEL
        if init:
            self.add_config_cmd(REG_MSLUTSEL, w[0] | (w[1] << 2) | (w[2] << 4) | (w[3] << 6)
                                | (x[0] << 8) | (x[1] << 16) | (x[2] << 24))
            if error:
                raise self.printer.config.error(error)
        else:
            self.set_register(REG_MSLUTSEL, w[0] | (w[1] << 2) | (w[2] << 4) | (w[3] << 6)
                              | (x[0] << 8) | (x[1] << 16) | (x[2] << 24))
            if error:
                logging.error(error)
                self.gcode.respond_info(error)
                return
            else:
                self.gcode.respond_info(success_msg)
        logging.info(success_msg)
    cmd_TMC_SET_WAVE_help = "Set wave correction factor for TMC2130 driver"
    def cmd_TMC_SET_WAVE(self, params):
        if 'FACTOR' in params:
            self.set_wave(self.gcode.get_float('FACTOR', params))

# Endstop wrapper that enables tmc2130 "sensorless homing"
class TMC2130VirtualEndstop:
    def __init__(self, tmc2130):
        self.tmc2130 = tmc2130
        if tmc2130.diag1_pin is None:
            raise pins.error("tmc2130 virtual endstop requires diag1_pin")
        ppins = tmc2130.printer.lookup_object('pins')
        self.mcu_endstop = ppins.setup_pin('endstop', tmc2130.diag1_pin)
        if self.mcu_endstop.get_mcu() is not tmc2130.mcu:
            raise pins.error("tmc2130 virtual endstop must be on same mcu")
        # Wrappers
        self.get_mcu = self.mcu_endstop.get_mcu
        self.add_stepper = self.mcu_endstop.add_stepper
        self.get_steppers = self.mcu_endstop.get_steppers
        self.home_start = self.mcu_endstop.home_start
        self.home_wait = self.mcu_endstop.home_wait
        self.query_endstop = self.mcu_endstop.query_endstop
        self.query_endstop_wait = self.mcu_endstop.query_endstop_wait
        self.TimeoutError = self.mcu_endstop.TimeoutError
    def home_prepare(self):
        gconf = self.tmc2130.reg_GCONF
        gconf &= ~GCONF_EN_PWM_MODE
        gconf |= GCONF_DIAG1_STALL
        self.tmc2130.set_register(REG_GCONF, gconf)
        self.tmc2130.set_register(REG_TCOOLTHRS, 0xfffff)
        self.mcu_endstop.home_prepare()
    def home_finalize(self):
        self.tmc2130.set_register(REG_GCONF, self.tmc2130.reg_GCONF)
        self.tmc2130.set_register(REG_TCOOLTHRS, 0)
        self.mcu_endstop.home_finalize()

def load_config_prefix(config):
    return TMC2130(config)
