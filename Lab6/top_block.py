#!/usr/bin/env python
##################################################
# Gnuradio Python Flow Graph
# Title: Top Block
# Generated: Wed Oct 26 01:27:55 2016
##################################################

from gnuradio import analog
from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import wxgui
from gnuradio.eng_option import eng_option
from gnuradio.fft import window
from gnuradio.filter import firdes
from gnuradio.wxgui import fftsink2
from gnuradio.wxgui import forms
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
        self.dn = dn = 10e-6
        self.dA = dA = 0.01

        ##################################################
        # Blocks
        ##################################################
        _dn_sizer = wx.BoxSizer(wx.VERTICAL)
        self._dn_text_box = forms.text_box(
        	parent=self.GetWin(),
        	sizer=_dn_sizer,
        	value=self.dn,
        	callback=self.set_dn,
        	label='dn',
        	converter=forms.float_converter(),
        	proportion=0,
        )
        self._dn_slider = forms.slider(
        	parent=self.GetWin(),
        	sizer=_dn_sizer,
        	value=self.dn,
        	callback=self.set_dn,
        	minimum=0,
        	maximum=30e-3,
        	num_steps=1000,
        	style=wx.SL_HORIZONTAL,
        	cast=float,
        	proportion=1,
        )
        self.Add(_dn_sizer)
        _dA_sizer = wx.BoxSizer(wx.VERTICAL)
        self._dA_text_box = forms.text_box(
        	parent=self.GetWin(),
        	sizer=_dA_sizer,
        	value=self.dA,
        	callback=self.set_dA,
        	label='dA',
        	converter=forms.float_converter(),
        	proportion=0,
        )
        self._dA_slider = forms.slider(
        	parent=self.GetWin(),
        	sizer=_dA_sizer,
        	value=self.dA,
        	callback=self.set_dA,
        	minimum=0.01,
        	maximum=1,
        	num_steps=100,
        	style=wx.SL_HORIZONTAL,
        	cast=float,
        	proportion=1,
        )
        self.Add(_dA_sizer)
        self.wxgui_fftsink2_0 = fftsink2.fft_sink_f(
        	self.GetWin(),
        	baseband_freq=0,
        	y_per_div=10,
        	y_divs=10,
        	ref_level=0,
        	ref_scale=2.0,
        	sample_rate=samp_rate,
        	fft_size=1024,
        	fft_rate=15,
        	average=False,
        	avg_alpha=None,
        	title="FFT Plot",
        	peak_hold=False,
        )
        self.Add(self.wxgui_fftsink2_0.win)
        self.blocks_transcendental_0 = blocks.transcendental("tanh", "float")
        self.blocks_throttle_0 = blocks.throttle(gr.sizeof_float*1, samp_rate)
        self.blocks_multiply_const_vxx_2 = blocks.multiply_const_vff((dn, ))
        self.blocks_multiply_const_vxx_0 = blocks.multiply_const_vff((20, ))
        self.blocks_file_source_0 = blocks.file_source(gr.sizeof_float*1, "/media/debarnab/Windows8_OS/Users/user/Desktop/5th Sem/EE340/Lab6/Input (1).bin", True)
        self.blocks_add_xx_0 = blocks.add_vff(1)
        self.analog_sig_source_x_0 = analog.sig_source_f(samp_rate, analog.GR_COS_WAVE, 3.3e3, dA, 0)

        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_transcendental_0, 0), (self.blocks_throttle_0, 0))
        self.connect((self.blocks_throttle_0, 0), (self.wxgui_fftsink2_0, 0))
        self.connect((self.blocks_file_source_0, 0), (self.blocks_multiply_const_vxx_2, 0))
        self.connect((self.blocks_multiply_const_vxx_0, 0), (self.blocks_transcendental_0, 0))
        self.connect((self.blocks_add_xx_0, 0), (self.blocks_multiply_const_vxx_0, 0))
        self.connect((self.blocks_multiply_const_vxx_2, 0), (self.blocks_add_xx_0, 0))
        self.connect((self.analog_sig_source_x_0, 0), (self.blocks_add_xx_0, 1))


# QT sink close method reimplementation

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.analog_sig_source_x_0.set_sampling_freq(self.samp_rate)
        self.wxgui_fftsink2_0.set_sample_rate(self.samp_rate)
        self.blocks_throttle_0.set_sample_rate(self.samp_rate)

    def get_dn(self):
        return self.dn

    def set_dn(self, dn):
        self.dn = dn
        self._dn_slider.set_value(self.dn)
        self._dn_text_box.set_value(self.dn)
        self.blocks_multiply_const_vxx_2.set_k((self.dn, ))

    def get_dA(self):
        return self.dA

    def set_dA(self, dA):
        self.dA = dA
        self._dA_slider.set_value(self.dA)
        self._dA_text_box.set_value(self.dA)
        self.analog_sig_source_x_0.set_amplitude(self.dA)

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

