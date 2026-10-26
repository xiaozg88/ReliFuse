# -*- coding: utf-8 -*-
"""Minimal OMML (Word equation) builders for the v31 methods equations."""

from xml.sax.saxutils import escape

M = 'xmlns:m="http://schemas.openxmlformats.org/officeDocument/2006/math"'


def r(t):
    return f'<m:r><m:t xml:space="preserve">{escape(t)}</m:t></m:r>'


def sub(base, sub):
    return f'<m:sSub><m:e>{base}</m:e><m:sub>{sub}</m:sub></m:sSub>'


def sup(base, sup_):
    return f'<m:sSup><m:e>{base}</m:e><m:sup>{sup_}</m:sup></m:sSup>'


def frac(num, den):
    return f'<m:f><m:num>{num}</m:num><m:den>{den}</m:den></m:f>'


def nary(body, lo, hi, char='∑'):
    return (f'<m:nary><m:naryPr><m:chr m:val="{char}"/></m:naryPr>'
            f'<m:sub>{lo}</m:sub><m:sup>{hi}</m:sup><m:e>{body}</m:e></m:nary>')


def delim(body, beg='[', end=']'):
    return (f'<m:d><m:dPr><m:begChr m:val="{escape(beg)}"/>'
            f'<m:endChr m:val="{escape(end)}"/></m:dPr><m:e>{body}</m:e></m:d>')


def eq(inner):
    return f'<m:oMath {M}>{inner}</m:oMath>'


EQ = {
    1: eq(sub(r('κ'), r('1')) + r('(i,j) = ') +
          frac(r('1 + cos(x_i, x_j)'), r('2'))),
    2: eq(sub(r('κ'), r('2')) + r('(i,j) = ') +
          frac(r('1'), r('1 + |i − j|'))),
    3: eq(sub(r('κ'), r('3')) + r('(i,j) = ') +
          frac(r('min(d_i, d_j)'), r('max(d_i, d_j)'))),
    4: eq(sub(r('κ'), r('4')) + r('(i,j) = exp(−‖x_i − x_j‖² / d)')),
    5: eq(sub(sup(r('g'), r('(l)')), r('e')) + r(' = β ') +
          nary(sub(r('α'), r('l,b,k')) + r(' κ') + sub(r('k'), r('(e)')),
               r('k=1'), r('4')) + r(' + (1 − β)')),
    6: eq(sub(r('z'), r('mc')) + r(' = ') + frac(sub(r('z'), r('m')),
                                                   sub(r('T'), r('m'))) +
          r(',  m ∈ {1, 2}')),
    7: eq(sub(r('p'), r('m')) + r(' = σ(') + sub(r('z'), r('mc')) + r(') = ') +
          frac(r('1'), r('1 + exp(−z_mc)'))),
    8: eq(r('H(p) = −[p ln p + (1 − p) ln(1 − p)]')),
    9: eq(
        r('agreement = ') + sub(r('p'), r('1')) + sub(r('p'), r('2'))
        + r(' + (1 − ') + sub(r('p'), r('1')) + r(')(1 − ') + sub(r('p'), r('2')) + r(')')
        + r(',  disagreement = |') + sub(r('p'), r('1')) + r(' − ') + sub(r('p'), r('2')) + r('|')),
    10: eq(sub(r('x'), r('gate')) + r(' = ') + delim(', '.join([
        r('z'), r('p'), r('H'), r('agreement'), r('disagreement'), sub(r('s'), r('struct')), r('ρ')
    ]))),
    11: eq(r('g = σ(') + sub(r('MLP'), r('g')) + r('(h)),  b = ') + sub(r('MLP'), r('b')) + r('(h)')),
    12: eq(sub(r('z'), r('fuse')) + r(' = ') + sub(r('z'), r('1c')) + r(' + g(') + sub(r('z'), r('2c')) + r(' − ') + sub(r('z'), r('1c')) + r(') + b')),
}