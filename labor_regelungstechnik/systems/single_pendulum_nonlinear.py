import control as ct
import numpy as np
from numpy import cos, sin, sqrt, tan


def system(t, states, inputs, params):
    # ~ unpacking the params
    m_x = params.get('m_x', 40)
    m_y = params.get('m_y', 3)
    y_max = params.get('y_max', 1.2)
    g = params.get('g', 9.81)
    c_varphi = params.get('c_varphi', 0.12)
    c_x = params.get('c_x', 0.5)
    k_x = params.get('k_x', 250)
    k_l = params.get('k_l', 500)
    l_0 = params.get('l_0', 0.23)

    k_vx = params.get('k_vx', 3.6)
    k_vl = params.get('k_vl', -1.65)

    # ~ unpacking the state
    L = states[0]
    X = states[1]
    l = states[2]
    phi = states[3]
    varphi = states[4]
    x = states[5]

    v_x = k_vx * inputs[0]
    v_l = k_vl * inputs[1]

    if x < 0 or x > 2.5:
        v_x = 0

    if l < 0 or l > 1.3:
        v_l = 0

    pass

    # ~ the main system equations
    d_L = L * k_l * l * m_x / (-l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
        varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
        varphi) ** 2 + l_0 * m_y ** 2 * cos(varphi) ** 2 - 2 * l_0 * m_y ** 2) - L * k_l * l * m_y * cos(
        varphi) ** 2 / (-l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
        varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
        varphi) ** 2 + l_0 * m_y ** 2 * cos(varphi) ** 2 - 2 * l_0 * m_y ** 2) + 2 * L * k_l * l * m_y / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) + L * k_l * l_0 * m_x / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) - L * k_l * l_0 * m_y * cos(varphi) ** 2 / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) + 2 * L * k_l * l_0 * m_y / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) + X * k_x * l * m_y * sin(varphi) / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) + X * k_x * l_0 * m_y * sin(varphi) / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) + c_varphi * m_y * phi * sin(varphi) * cos(varphi) / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) + g * l * m_y ** 2 * sin(varphi) ** 2 * cos(varphi) / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) + g * l_0 * m_y ** 2 * sin(varphi) ** 2 * cos(
        varphi) / (-l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
        varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
        varphi) ** 2 + l_0 * m_y ** 2 * cos(varphi) ** 2 - 2 * l_0 * m_y ** 2) - k_l * l * m_x * v_l / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) + k_l * l * m_y * v_l * cos(varphi) ** 2 / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) - 2 * k_l * l * m_y * v_l / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) - k_l * l_0 * m_x * v_l / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) + k_l * l_0 * m_y * v_l * cos(varphi) ** 2 / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) - 2 * k_l * l_0 * m_y * v_l / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) - k_x * l * m_y * v_x * sin(varphi) / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) - k_x * l_0 * m_y * v_x * sin(varphi) / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) - l ** 2 * m_x * m_y * phi ** 2 / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) + l ** 2 * m_y ** 2 * phi ** 2 * sin(varphi) ** 2 / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) + l ** 2 * m_y ** 2 * phi ** 2 * cos(varphi) ** 2 / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) - 2 * l ** 2 * m_y ** 2 * phi ** 2 / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) - 2 * l * l_0 * m_x * m_y * phi ** 2 / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) + 2 * l * l_0 * m_y ** 2 * phi ** 2 * sin(
        varphi) ** 2 / (-l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
        varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
        varphi) ** 2 + l_0 * m_y ** 2 * cos(
        varphi) ** 2 - 2 * l_0 * m_y ** 2) + 2 * l * l_0 * m_y ** 2 * phi ** 2 * cos(varphi) ** 2 / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) - 4 * l * l_0 * m_y ** 2 * phi ** 2 / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) - l_0 ** 2 * m_x * m_y * phi ** 2 / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) + l_0 ** 2 * m_y ** 2 * phi ** 2 * sin(varphi) ** 2 / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) + l_0 ** 2 * m_y ** 2 * phi ** 2 * cos(varphi) ** 2 / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l_0 * m_y ** 2) - 2 * l_0 ** 2 * m_y ** 2 * phi ** 2 / (
                      -l * m_x * m_y + l * m_y ** 2 * sin(varphi) ** 2 + l * m_y ** 2 * cos(
                  varphi) ** 2 - 2 * l * m_y ** 2 - l_0 * m_x * m_y + l_0 * m_y ** 2 * sin(
                  varphi) ** 2 + l_0 * m_y ** 2 * cos(varphi) ** 2 - 2 * l_0 * m_y ** 2)
    d_X = L * k_l * l * sin(varphi) / (-l * m_x + l * m_y * sin(varphi) ** 2 + l * m_y * cos(
        varphi) ** 2 - 2 * l * m_y - l_0 * m_x + l_0 * m_y * sin(varphi) ** 2 + l_0 * m_y * cos(
        varphi) ** 2 - 2 * l_0 * m_y) + L * k_l * l_0 * sin(varphi) / (
                      -l * m_x + l * m_y * sin(varphi) ** 2 + l * m_y * cos(
                  varphi) ** 2 - 2 * l * m_y - l_0 * m_x + l_0 * m_y * sin(varphi) ** 2 + l_0 * m_y * cos(
                  varphi) ** 2 - 2 * l_0 * m_y) + X * k_x * l / (
                      -l * m_x + l * m_y * sin(varphi) ** 2 + l * m_y * cos(
                  varphi) ** 2 - 2 * l * m_y - l_0 * m_x + l_0 * m_y * sin(varphi) ** 2 + l_0 * m_y * cos(
                  varphi) ** 2 - 2 * l_0 * m_y) + X * k_x * l_0 / (
                      -l * m_x + l * m_y * sin(varphi) ** 2 + l * m_y * cos(
                  varphi) ** 2 - 2 * l * m_y - l_0 * m_x + l_0 * m_y * sin(varphi) ** 2 + l_0 * m_y * cos(
                  varphi) ** 2 - 2 * l_0 * m_y) + c_varphi * phi * cos(varphi) / (
                      -l * m_x + l * m_y * sin(varphi) ** 2 + l * m_y * cos(
                  varphi) ** 2 - 2 * l * m_y - l_0 * m_x + l_0 * m_y * sin(varphi) ** 2 + l_0 * m_y * cos(
                  varphi) ** 2 - 2 * l_0 * m_y) + g * l * m_y * sin(varphi) * cos(varphi) / (
                      -l * m_x + l * m_y * sin(varphi) ** 2 + l * m_y * cos(
                  varphi) ** 2 - 2 * l * m_y - l_0 * m_x + l_0 * m_y * sin(varphi) ** 2 + l_0 * m_y * cos(
                  varphi) ** 2 - 2 * l_0 * m_y) + g * l_0 * m_y * sin(varphi) * cos(varphi) / (
                      -l * m_x + l * m_y * sin(varphi) ** 2 + l * m_y * cos(
                  varphi) ** 2 - 2 * l * m_y - l_0 * m_x + l_0 * m_y * sin(varphi) ** 2 + l_0 * m_y * cos(
                  varphi) ** 2 - 2 * l_0 * m_y) - k_l * l * v_l * sin(varphi) / (
                      -l * m_x + l * m_y * sin(varphi) ** 2 + l * m_y * cos(
                  varphi) ** 2 - 2 * l * m_y - l_0 * m_x + l_0 * m_y * sin(varphi) ** 2 + l_0 * m_y * cos(
                  varphi) ** 2 - 2 * l_0 * m_y) - k_l * l_0 * v_l * sin(varphi) / (
                      -l * m_x + l * m_y * sin(varphi) ** 2 + l * m_y * cos(
                  varphi) ** 2 - 2 * l * m_y - l_0 * m_x + l_0 * m_y * sin(varphi) ** 2 + l_0 * m_y * cos(
                  varphi) ** 2 - 2 * l_0 * m_y) - k_x * l * v_x / (
                      -l * m_x + l * m_y * sin(varphi) ** 2 + l * m_y * cos(
                  varphi) ** 2 - 2 * l * m_y - l_0 * m_x + l_0 * m_y * sin(varphi) ** 2 + l_0 * m_y * cos(
                  varphi) ** 2 - 2 * l_0 * m_y) - k_x * l_0 * v_x / (
                      -l * m_x + l * m_y * sin(varphi) ** 2 + l * m_y * cos(
                  varphi) ** 2 - 2 * l * m_y - l_0 * m_x + l_0 * m_y * sin(varphi) ** 2 + l_0 * m_y * cos(
                  varphi) ** 2 - 2 * l_0 * m_y)
    d_l = L
    d_phi = L * k_l * l * m_y * sin(varphi) * cos(varphi) / (
                -l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
            varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
            varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
            varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
            varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
            varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) + L * k_l * l_0 * m_y * sin(varphi) * cos(varphi) / (
                        -l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
                    varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
                    varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
                    varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) + 2 * L * l * m_x * m_y * phi / (
                        -l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
                    varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
                    varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
                    varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) - 2 * L * l * m_y ** 2 * phi * sin(
        varphi) ** 2 / (-l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
        varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
        varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
        varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
        varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
        varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) - 2 * L * l * m_y ** 2 * phi * cos(varphi) ** 2 / (
                        -l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
                    varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
                    varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
                    varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) + 4 * L * l * m_y ** 2 * phi / (
                        -l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
                    varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
                    varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
                    varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) + 2 * L * l_0 * m_x * m_y * phi / (
                        -l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
                    varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
                    varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
                    varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) - 2 * L * l_0 * m_y ** 2 * phi * sin(
        varphi) ** 2 / (-l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
        varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
        varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
        varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
        varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
        varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) - 2 * L * l_0 * m_y ** 2 * phi * cos(varphi) ** 2 / (
                        -l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
                    varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
                    varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
                    varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) + 4 * L * l_0 * m_y ** 2 * phi / (
                        -l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
                    varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
                    varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
                    varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) + X * k_x * l * m_y * cos(varphi) / (
                        -l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
                    varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
                    varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
                    varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) + X * k_x * l_0 * m_y * cos(varphi) / (
                        -l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
                    varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
                    varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
                    varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) + c_varphi * m_x * phi / (
                        -l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
                    varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
                    varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
                    varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) - c_varphi * m_y * phi * sin(varphi) ** 2 / (
                        -l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
                    varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
                    varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
                    varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) + 2 * c_varphi * m_y * phi / (
                        -l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
                    varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
                    varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
                    varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) + g * l * m_x * m_y * sin(varphi) / (
                        -l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
                    varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
                    varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
                    varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) - g * l * m_y ** 2 * sin(varphi) ** 3 / (
                        -l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
                    varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
                    varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
                    varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) + 2 * g * l * m_y ** 2 * sin(varphi) / (
                        -l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
                    varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
                    varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
                    varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) + g * l_0 * m_x * m_y * sin(varphi) / (
                        -l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
                    varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
                    varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
                    varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) - g * l_0 * m_y ** 2 * sin(varphi) ** 3 / (
                        -l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
                    varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
                    varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
                    varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) + 2 * g * l_0 * m_y ** 2 * sin(varphi) / (
                        -l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
                    varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
                    varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
                    varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) - k_l * l * m_y * v_l * sin(varphi) * cos(
        varphi) / (-l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
        varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
        varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
        varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
        varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
        varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) - k_l * l_0 * m_y * v_l * sin(varphi) * cos(varphi) / (
                        -l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
                    varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
                    varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
                    varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) - k_x * l * m_y * v_x * cos(varphi) / (
                        -l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
                    varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
                    varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
                    varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2) - k_x * l_0 * m_y * v_x * cos(varphi) / (
                        -l ** 2 * m_x * m_y + l ** 2 * m_y ** 2 * sin(varphi) ** 2 + l ** 2 * m_y ** 2 * cos(
                    varphi) ** 2 - 2 * l ** 2 * m_y ** 2 - 2 * l * l_0 * m_x * m_y + 2 * l * l_0 * m_y ** 2 * sin(
                    varphi) ** 2 + 2 * l * l_0 * m_y ** 2 * cos(
                    varphi) ** 2 - 4 * l * l_0 * m_y ** 2 - l_0 ** 2 * m_x * m_y + l_0 ** 2 * m_y ** 2 * sin(
                    varphi) ** 2 + l_0 ** 2 * m_y ** 2 * cos(varphi) ** 2 - 2 * l_0 ** 2 * m_y ** 2)
    d_varphi = phi
    d_x = X

    return [d_L, d_X, d_l, d_phi, d_varphi, d_x]


def output(t, states, inputs, params):
    m_x = params.get('m_x', 30)
    m_y = params.get('m_y', 3)
    y_max = params.get('y_max', 1.2)
    g = params.get('g', 9.81)
    c_varphi = params.get('c_varphi', 0.5)
    c_x = params.get('c_x', 0.5)
    k_x = params.get('k_x', 200)
    k_l = params.get('k_l', 200)
    l_0 = params.get('l_0', 0.1)
    k_vx = params.get('k_vx', 4)
    k_vl = params.get('k_vl', 2)
    k_phi = params.get('k_phi', 0.5)

    # ~ unpacking the state
    L = states[0]
    X = states[1]
    l = states[2]
    phi = states[3]
    varphi = states[4]
    x = states[5]

    return [x, l, k_phi * np.degrees(varphi)]


output_names = ('x', 'l', 'varphi')
io_system = ct.NonlinearIOSystem(
    system, output,
    inputs=('v_x', 'v_l'),
    outputs=output_names,
    states=('L', 'X', 'l', 'phi', 'varphi', 'x'),
    name='system',
)