#!/usr/bin/env python
##################################################
# Gnuradio Python Flow Graph
# Title: Top Block
# Generated: Wed Oct 26 13:07:47 2016
##################################################

from gnuradio import analog
from gnuradio import blocks
from gnuradio import digital
from gnuradio import eng_notation
from gnuradio import filter
from gnuradio import gr
from gnuradio import wxgui
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from gnuradio.filter import pfb
from gnuradio.wxgui import forms
from gnuradio.wxgui import scopesink2
from grc_gnuradio import wxgui as grc_wxgui
from math import pi
from optparse import OptionParser
import numpy
import wx

class top_block(grc_wxgui.top_block_gui):

    def __init__(self):
        grc_wxgui.top_block_gui.__init__(self, title="Top Block")
        _icon_path = "/usr/share/icons/hicolor/32x32/apps/gnuradio-grc.png"
        self.SetIcon(wx.Icon(_icon_path, wx.BITMAP_TYPE_ANY))

        ##################################################
        # Variables
        ##################################################
        self.sps = sps = 4
        self.nfilts = nfilts = 32
        self.ntaps = ntaps = 11*nfilts*sps
        self.excess_bw = excess_bw = 400e-3
        self.tx_taps = tx_taps = firdes.root_raised_cosine(nfilts,nfilts,1.0,excess_bw,ntaps)
        self.timing_bw = timing_bw = 2*pi/100
        self.tap4 = tap4 = 0.0625
        self.tap3 = tap3 = -0.125
        self.tap2 = tap2 = 0.25
        self.tap1 = tap1 = -0.5
        self.sqrt_2 = sqrt_2 = 1.41421356237
        self.samp_rate = samp_rate = 32000
        self.rx_taps = rx_taps = filter.firdes.root_raised_cosine(nfilts, nfilts*sps, 1.0,
        excess_bw, ntaps)
        self.const_points = const_points = 8

        ##################################################
        # Blocks
        ##################################################
        _tap4_sizer = wx.BoxSizer(wx.VERTICAL)
        self._tap4_text_box = forms.text_box(
        	parent=self.GetWin(),
        	sizer=_tap4_sizer,
        	value=self.tap4,
        	callback=self.set_tap4,
        	label='tap4',
        	converter=forms.float_converter(),
        	proportion=0,
        )
        self._tap4_slider = forms.slider(
        	parent=self.GetWin(),
        	sizer=_tap4_sizer,
        	value=self.tap4,
        	callback=self.set_tap4,
        	minimum=-1,
        	maximum=1,
        	num_steps=1000,
        	style=wx.SL_HORIZONTAL,
        	cast=float,
        	proportion=1,
        )
        self.Add(_tap4_sizer)
        _tap3_sizer = wx.BoxSizer(wx.VERTICAL)
        self._tap3_text_box = forms.text_box(
        	parent=self.GetWin(),
        	sizer=_tap3_sizer,
        	value=self.tap3,
        	callback=self.set_tap3,
        	label='tap3',
        	converter=forms.float_converter(),
        	proportion=0,
        )
        self._tap3_slider = forms.slider(
        	parent=self.GetWin(),
        	sizer=_tap3_sizer,
        	value=self.tap3,
        	callback=self.set_tap3,
        	minimum=-1,
        	maximum=1,
        	num_steps=1000,
        	style=wx.SL_HORIZONTAL,
        	cast=float,
        	proportion=1,
        )
        self.Add(_tap3_sizer)
        _tap2_sizer = wx.BoxSizer(wx.VERTICAL)
        self._tap2_text_box = forms.text_box(
        	parent=self.GetWin(),
        	sizer=_tap2_sizer,
        	value=self.tap2,
        	callback=self.set_tap2,
        	label='tap2',
        	converter=forms.float_converter(),
        	proportion=0,
        )
        self._tap2_slider = forms.slider(
        	parent=self.GetWin(),
        	sizer=_tap2_sizer,
        	value=self.tap2,
        	callback=self.set_tap2,
        	minimum=-1,
        	maximum=1,
        	num_steps=1000,
        	style=wx.SL_HORIZONTAL,
        	cast=float,
        	proportion=1,
        )
        self.Add(_tap2_sizer)
        _tap1_sizer = wx.BoxSizer(wx.VERTICAL)
        self._tap1_text_box = forms.text_box(
        	parent=self.GetWin(),
        	sizer=_tap1_sizer,
        	value=self.tap1,
        	callback=self.set_tap1,
        	label='tap1',
        	converter=forms.float_converter(),
        	proportion=0,
        )
        self._tap1_slider = forms.slider(
        	parent=self.GetWin(),
        	sizer=_tap1_sizer,
        	value=self.tap1,
        	callback=self.set_tap1,
        	minimum=-1,
        	maximum=1,
        	num_steps=1000,
        	style=wx.SL_HORIZONTAL,
        	cast=float,
        	proportion=1,
        )
        self.Add(_tap1_sizer)
        self.notebook_0 = self.notebook_0 = wx.Notebook(self.GetWin(), style=wx.NB_TOP)
        self.notebook_0.AddPage(grc_wxgui.Panel(self.notebook_0), "tab1")
        self.notebook_0.AddPage(grc_wxgui.Panel(self.notebook_0), "tab2")
        self.notebook_0.AddPage(grc_wxgui.Panel(self.notebook_0), "tab3")
        self.notebook_0.AddPage(grc_wxgui.Panel(self.notebook_0), "tab4")
        self.Add(self.notebook_0)
        self.wxgui_scopesink2_0 = scopesink2.scope_sink_c(
        	self.notebook_0.GetPage(0).GetWin(),
        	title="Scope Plot",
        	sample_rate=samp_rate,
        	v_scale=0,
        	v_offset=0,
        	t_scale=0,
        	ac_couple=False,
        	xy_mode=False,
        	num_inputs=1,
        	trig_mode=wxgui.TRIG_MODE_AUTO,
        	y_axis_label="Counts",
        )
        self.notebook_0.GetPage(0).Add(self.wxgui_scopesink2_0.win)
        self.pfb_arb_resampler_xxx_0 = pfb.arb_resampler_ccf(
        	  sps,
                  taps=(tx_taps),
        	  flt_size=32)
        	
        self.iir_filter_ffd_0_0 = filter.iir_filter_ffd(([1,tap1,tap2,tap3,tap4]), ([1]), True)
        self.iir_filter_ffd_0 = filter.iir_filter_ffd(([1,tap1,tap2,tap3,tap4]), ([1]), True)
        self.digital_pfb_clock_sync_xxx_0 = digital.pfb_clock_sync_ccf(sps, timing_bw, (rx_taps), 32, 16, 1.5, 1)
        self.digital_chunks_to_symbols_xx_0 = digital.chunks_to_symbols_bc((((1),(0.707+0.707j),(1j),(-0.707+0.707j),(-1),(-0.707-0.707j),(-1j),(0.707-0.707j))), 1)
        self.blocks_throttle_0 = blocks.throttle(gr.sizeof_gr_complex*1, samp_rate*sps)
        self.blocks_multiply_const_vxx_0 = blocks.multiply_const_vcc((0.5, ))
        self.blocks_float_to_complex_0 = blocks.float_to_complex(1)
        self.blocks_delay_0 = blocks.delay(gr.sizeof_gr_complex*1, 1)
        self.blocks_complex_to_float_0 = blocks.complex_to_float(1)
        self.blocks_add_xx_0 = blocks.add_vcc(1)
        self.analog_random_source_x_0 = blocks.vector_source_b(map(int, numpy.random.randint(0, const_points, 1000)), True)
        self.analog_noise_source_x_0 = analog.noise_source_c(analog.GR_GAUSSIAN, 0.07, 0)

        ##################################################
        # Connections
        ##################################################
        self.connect((self.analog_random_source_x_0, 0), (self.digital_chunks_to_symbols_xx_0, 0))
        self.connect((self.blocks_multiply_const_vxx_0, 0), (self.blocks_delay_0, 0))
        self.connect((self.blocks_add_xx_0, 0), (self.blocks_complex_to_float_0, 0))
        self.connect((self.blocks_complex_to_float_0, 1), (self.iir_filter_ffd_0_0, 0))
        self.connect((self.blocks_complex_to_float_0, 0), (self.iir_filter_ffd_0, 0))
        self.connect((self.iir_filter_ffd_0, 0), (self.blocks_float_to_complex_0, 0))
        self.connect((self.iir_filter_ffd_0_0, 0), (self.blocks_float_to_complex_0, 1))
        self.connect((self.blocks_float_to_complex_0, 0), (self.pfb_arb_resampler_xxx_0, 0))
        self.connect((self.blocks_throttle_0, 0), (self.digital_pfb_clock_sync_xxx_0, 0))
        self.connect((self.digital_pfb_clock_sync_xxx_0, 0), (self.wxgui_scopesink2_0, 0))
        self.connect((self.pfb_arb_resampler_xxx_0, 0), (self.blocks_throttle_0, 0))
        self.connect((self.digital_chunks_to_symbols_xx_0, 0), (self.blocks_add_xx_0, 0))
        self.connect((self.blocks_delay_0, 0), (self.blocks_add_xx_0, 1))
        self.connect((self.digital_chunks_to_symbols_xx_0, 0), (self.blocks_multiply_const_vxx_0, 0))
        self.connect((self.analog_noise_source_x_0, 0), (self.blocks_add_xx_0, 2))


# QT sink close method reimplementation

    def get_sps(self):
        return self.sps

    def set_sps(self, sps):
        self.sps = sps
        self.set_rx_taps(filter.firdes.root_raised_cosine(self.nfilts, self.nfilts*self.sps, 1.0,
        self.excess_bw, self.ntaps))
        self.set_ntaps(11*self.nfilts*self.sps)
        self.pfb_arb_resampler_xxx_0.set_rate(self.sps)
        self.blocks_throttle_0.set_sample_rate(self.samp_rate*self.sps)

    def get_nfilts(self):
        return self.nfilts

    def set_nfilts(self, nfilts):
        self.nfilts = nfilts
        self.set_rx_taps(filter.firdes.root_raised_cosine(self.nfilts, self.nfilts*self.sps, 1.0,
        self.excess_bw, self.ntaps))
        self.set_tx_taps(firdes.root_raised_cosine(self.nfilts,self.nfilts,1.0,self.excess_bw,self.ntaps))
        self.set_ntaps(11*self.nfilts*self.sps)

    def get_ntaps(self):
        return self.ntaps

    def set_ntaps(self, ntaps):
        self.ntaps = ntaps
        self.set_rx_taps(filter.firdes.root_raised_cosine(self.nfilts, self.nfilts*self.sps, 1.0,
        self.excess_bw, self.ntaps))
        self.set_tx_taps(firdes.root_raised_cosine(self.nfilts,self.nfilts,1.0,self.excess_bw,self.ntaps))

    def get_excess_bw(self):
        return self.excess_bw

    def set_excess_bw(self, excess_bw):
        self.excess_bw = excess_bw
        self.set_rx_taps(filter.firdes.root_raised_cosine(self.nfilts, self.nfilts*self.sps, 1.0,
        self.excess_bw, self.ntaps))
        self.set_tx_taps(firdes.root_raised_cosine(self.nfilts,self.nfilts,1.0,self.excess_bw,self.ntaps))

    def get_tx_taps(self):
        return self.tx_taps

    def set_tx_taps(self, tx_taps):
        self.tx_taps = tx_taps
        self.pfb_arb_resampler_xxx_0.set_taps((self.tx_taps))

    def get_timing_bw(self):
        return self.timing_bw

    def set_timing_bw(self, timing_bw):
        self.timing_bw = timing_bw
        self.digital_pfb_clock_sync_xxx_0.set_loop_bandwidth(self.timing_bw)

    def get_tap4(self):
        return self.tap4

    def set_tap4(self, tap4):
        self.tap4 = tap4
        self.iir_filter_ffd_0_0.set_taps(([1,self.tap1,self.tap2,self.tap3,self.tap4]), ([1]))
        self._tap4_slider.set_value(self.tap4)
        self._tap4_text_box.set_value(self.tap4)
        self.iir_filter_ffd_0.set_taps(([1,self.tap1,self.tap2,self.tap3,self.tap4]), ([1]))

    def get_tap3(self):
        return self.tap3

    def set_tap3(self, tap3):
        self.tap3 = tap3
        self.iir_filter_ffd_0_0.set_taps(([1,self.tap1,self.tap2,self.tap3,self.tap4]), ([1]))
        self._tap3_slider.set_value(self.tap3)
        self._tap3_text_box.set_value(self.tap3)
        self.iir_filter_ffd_0.set_taps(([1,self.tap1,self.tap2,self.tap3,self.tap4]), ([1]))

    def get_tap2(self):
        return self.tap2

    def set_tap2(self, tap2):
        self.tap2 = tap2
        self.iir_filter_ffd_0_0.set_taps(([1,self.tap1,self.tap2,self.tap3,self.tap4]), ([1]))
        self._tap2_slider.set_value(self.tap2)
        self._tap2_text_box.set_value(self.tap2)
        self.iir_filter_ffd_0.set_taps(([1,self.tap1,self.tap2,self.tap3,self.tap4]), ([1]))

    def get_tap1(self):
        return self.tap1

    def set_tap1(self, tap1):
        self.tap1 = tap1
        self.iir_filter_ffd_0_0.set_taps(([1,self.tap1,self.tap2,self.tap3,self.tap4]), ([1]))
        self._tap1_slider.set_value(self.tap1)
        self._tap1_text_box.set_value(self.tap1)
        self.iir_filter_ffd_0.set_taps(([1,self.tap1,self.tap2,self.tap3,self.tap4]), ([1]))

    def get_sqrt_2(self):
        return self.sqrt_2

    def set_sqrt_2(self, sqrt_2):
        self.sqrt_2 = sqrt_2

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.wxgui_scopesink2_0.set_sample_rate(self.samp_rate)
        self.blocks_throttle_0.set_sample_rate(self.samp_rate*self.sps)

    def get_rx_taps(self):
        return self.rx_taps

    def set_rx_taps(self, rx_taps):
        self.rx_taps = rx_taps
        self.digital_pfb_clock_sync_xxx_0.set_taps((self.rx_taps))

    def get_const_points(self):
        return self.const_points

    def set_const_points(self, const_points):
        self.const_points = const_points

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

