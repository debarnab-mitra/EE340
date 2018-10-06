#!/usr/bin/env python
##################################################
# Gnuradio Python Flow Graph
# Title: Top Block
# Generated: Wed Oct 26 00:57:26 2016
##################################################

from gnuradio import analog
from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import filter
from gnuradio import gr
from gnuradio import wxgui
from gnuradio.eng_option import eng_option
from gnuradio.fft import window
from gnuradio.filter import firdes
from gnuradio.wxgui import fftsink2
from gnuradio.wxgui import scopesink2
from grc_gnuradio import wxgui as grc_wxgui
from optparse import OptionParser
import wx

class top_block(grc_wxgui.top_block_gui):

    def __init__(self):
        grc_wxgui.top_block_gui.__init__(self, title="Top Block")
        _icon_path = "/usr/share/icons/hicolor/32x32/apps/gnuradio-grc.png"
        self.SetIcon(wx.Icon(_icon_path, wx.BITMAP_TYPE_ANY))

        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate = 32000

        ##################################################
        # Blocks
        ##################################################
        self.debu = self.debu = wx.Notebook(self.GetWin(), style=wx.NB_TOP)
        self.debu.AddPage(grc_wxgui.Panel(self.debu), "tab1")
        self.debu.AddPage(grc_wxgui.Panel(self.debu), "tab2")
        self.debu.AddPage(grc_wxgui.Panel(self.debu), "demod")
        self.debu.AddPage(grc_wxgui.Panel(self.debu), "demodft")
        self.Add(self.debu)
        self.wxgui_scopesink2_0_0 = scopesink2.scope_sink_f(
        	self.debu.GetPage(2).GetWin(),
        	title="Scope Plot",
        	sample_rate=40e3,
        	v_scale=0,
        	v_offset=0,
        	t_scale=0,
        	ac_couple=False,
        	xy_mode=False,
        	num_inputs=1,
        	trig_mode=wxgui.TRIG_MODE_AUTO,
        	y_axis_label="Counts",
        )
        self.debu.GetPage(2).Add(self.wxgui_scopesink2_0_0.win)
        self.wxgui_scopesink2_0 = scopesink2.scope_sink_f(
        	self.debu.GetPage(0).GetWin(),
        	title="Scope Plot",
        	sample_rate=640e3,
        	v_scale=0,
        	v_offset=0,
        	t_scale=0,
        	ac_couple=False,
        	xy_mode=False,
        	num_inputs=1,
        	trig_mode=wxgui.TRIG_MODE_AUTO,
        	y_axis_label="Counts",
        )
        self.debu.GetPage(0).Add(self.wxgui_scopesink2_0.win)
        self.wxgui_fftsink2_0_0 = fftsink2.fft_sink_f(
        	self.debu.GetPage(3).GetWin(),
        	baseband_freq=0,
        	y_per_div=10,
        	y_divs=10,
        	ref_level=0,
        	ref_scale=2.0,
        	sample_rate=40e3,
        	fft_size=1024,
        	fft_rate=15,
        	average=False,
        	avg_alpha=None,
        	title="FFT Plot",
        	peak_hold=False,
        )
        self.debu.GetPage(3).Add(self.wxgui_fftsink2_0_0.win)
        self.wxgui_fftsink2_0 = fftsink2.fft_sink_f(
        	self.debu.GetPage(1).GetWin(),
        	baseband_freq=0,
        	y_per_div=10,
        	y_divs=10,
        	ref_level=0,
        	ref_scale=2.0,
        	sample_rate=640e3,
        	fft_size=1024,
        	fft_rate=15,
        	average=False,
        	avg_alpha=None,
        	title="FFT Plot",
        	peak_hold=False,
        )
        self.debu.GetPage(1).Add(self.wxgui_fftsink2_0.win)
        self.rational_resampler_xxx_0 = filter.rational_resampler_fff(
                interpolation=16,
                decimation=1,
                taps=None,
                fractional_bw=None,
        )
        self.low_pass_filter_0 = filter.fir_filter_fff(16, firdes.low_pass(
        	1, 640e3, 20e3, 4000, firdes.WIN_HAMMING, 6.76))
        self.iir_filter_ffd_1_0 = filter.iir_filter_ffd((1, ), (1,0.95), True)
        self.iir_filter_ffd_1 = filter.iir_filter_ffd((1,-0.95), (1, ), True)
        self.iir_filter_ffd_0 = filter.iir_filter_ffd((25e-6, ), (1,1), True)
        self.blocks_throttle_0_0 = blocks.throttle(gr.sizeof_float*1, 40e3)
        self.blocks_throttle_0 = blocks.throttle(gr.sizeof_float*1, 40e3)
        self.blocks_multiply_xx_0 = blocks.multiply_vcc(1)
        self.blocks_multiply_conjugate_cc_0 = blocks.multiply_conjugate_cc(1)
        self.blocks_delay_0 = blocks.delay(gr.sizeof_gr_complex*1, 1)
        self.blocks_complex_to_real_0 = blocks.complex_to_real(1)
        self.blocks_complex_to_arg_0 = blocks.complex_to_arg(1)
        self.blocks_add_xx_1 = blocks.add_vcc(1)
        self.blocks_add_xx_0 = blocks.add_vff(1)
        self.analog_sig_source_x_2 = analog.sig_source_f(40e3, analog.GR_COS_WAVE, 1100, 0.5, 0)
        self.analog_sig_source_x_1 = analog.sig_source_c(640e3, analog.GR_COS_WAVE, 100000, 1, 0)
        self.analog_sig_source_x_0 = analog.sig_source_f(40e3, analog.GR_COS_WAVE, 11000, 0.5, 0)
        self.analog_phase_modulator_fc_0 = analog.phase_modulator_fc(471.24e3)
        self.analog_noise_source_x_0 = analog.noise_source_c(analog.GR_GAUSSIAN, 0.2, 0)

        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_multiply_conjugate_cc_0, 0), (self.blocks_complex_to_arg_0, 0))
        self.connect((self.blocks_add_xx_1, 0), (self.blocks_multiply_conjugate_cc_0, 0))
        self.connect((self.blocks_delay_0, 0), (self.blocks_multiply_conjugate_cc_0, 1))
        self.connect((self.blocks_add_xx_1, 0), (self.blocks_delay_0, 0))
        self.connect((self.analog_noise_source_x_0, 0), (self.blocks_add_xx_1, 1))
        self.connect((self.analog_sig_source_x_1, 0), (self.blocks_multiply_xx_0, 1))
        self.connect((self.blocks_complex_to_real_0, 0), (self.wxgui_fftsink2_0, 0))
        self.connect((self.blocks_complex_to_real_0, 0), (self.wxgui_scopesink2_0, 0))
        self.connect((self.blocks_add_xx_1, 0), (self.blocks_complex_to_real_0, 0))
        self.connect((self.blocks_multiply_xx_0, 0), (self.blocks_add_xx_1, 0))
        self.connect((self.analog_phase_modulator_fc_0, 0), (self.blocks_multiply_xx_0, 0))
        self.connect((self.blocks_add_xx_0, 0), (self.iir_filter_ffd_0, 0))
        self.connect((self.rational_resampler_xxx_0, 0), (self.analog_phase_modulator_fc_0, 0))
        self.connect((self.analog_sig_source_x_0, 0), (self.blocks_throttle_0_0, 0))
        self.connect((self.blocks_throttle_0_0, 0), (self.blocks_add_xx_0, 1))
        self.connect((self.analog_sig_source_x_2, 0), (self.blocks_throttle_0, 0))
        self.connect((self.blocks_throttle_0, 0), (self.blocks_add_xx_0, 0))
        self.connect((self.iir_filter_ffd_0, 0), (self.iir_filter_ffd_1, 0))
        self.connect((self.iir_filter_ffd_1, 0), (self.rational_resampler_xxx_0, 0))
        self.connect((self.blocks_complex_to_arg_0, 0), (self.low_pass_filter_0, 0))
        self.connect((self.low_pass_filter_0, 0), (self.iir_filter_ffd_1_0, 0))
        self.connect((self.iir_filter_ffd_1_0, 0), (self.wxgui_scopesink2_0_0, 0))
        self.connect((self.iir_filter_ffd_1_0, 0), (self.wxgui_fftsink2_0_0, 0))


# QT sink close method reimplementation

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate

if __name__ == '__main__':
    import ctypes
    import os
    if os.name == 'posix':
        try:
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except:
            print "Warning: failed to XInitThreads()"
    parser = OptionParser(option_class=eng_option, usage="%prog: [options]")
    (options, args) = parser.parse_args()
    tb = top_block()
    tb.Start(True)
    tb.Wait()

