import numpy as np
from numpy import sin, cos, arctan2
from itertools import cycle
from sys import argv, exit
import pyqtgraph as pg
from pyqtgraph import QtCore, QtWidgets, QtGui

MIN_DX = -15
MAX_DX = 15
MAX_THETA = np.pi
MIN_THETA = -np.pi
MAX_DTHETA = 1.1
MIN_DTHETA = 0

class InvertedPendulum(QtWidgets.QWidget):
    '''Inicjalizacja stałych:
    M - masa wózka
    m - masa kulki
    l - długość ramienia wahadła

    Warunków początkowych:
    x0 - początkowe położenie wózka
    dx0 - początkowa prędkość wózka
    theta0 - początkowe położenie wahadla
    dtheta0 - początkowa prędkość wahadła

    Zakłócenia zewnętrznego:
    dis_cyc - zmienna odpowiada za to, czy zakłócenie jest zapętlone
    disruption - wartości zakłócenia w kolejnych chwilach czasowych

    Parametry planszy/obrazka:
    iw, ih - szerokość i wysokość obrazka
    x_max - maksymalna współrzędna pozioma (oś x jest symetryczna, więc minimalna wynosi -x_max)
    h_min - minialna współrzędna pionowa
    h_max - maksymalna współrzędna pionowa

    Powyższe dane są pobierane z pliku jeśli zmienna f_name nie jest pusta'''
    def __init__(self, M=10, m=5, l=50, x0=0, theta0=0, dx0=0, dtheta0=0, dis_cyc=True, disruption=[0], iw=1000, ih=500, x_max=100, h_min=0, h_max=100, f_name=None):
        if f_name:
            with open(f_name) as f_handle:
                lines = f_handle.readlines()
                init_cond = lines[0].split(' ')
                self.M, self.m, self.l, self.x0, self.theta0, self.dx0, self.dtheta0 = [float(el) for el in init_cond[:7]]
                self.image_w, self.image_h, self.x_max, self.h_min, self.h_max = [int(el) for el in init_cond[-5:]]
                if lines[1].strip()=='1':
                    self.disruption = cycle([float(el) for el in lines[2].split(' ')])
                else:
                    self.disruption = iter([float(el) for el in lines[2].split(' ')])
        else:
            self.M, self.m, self.l, self.x0, self.theta0, self.dx0, self.dtheta0 = M, m, l, x0, theta0, dx0, dtheta0
            self.image_w, self.image_h, self.x_max, self.h_min, self.h_max = iw, ih, x_max, h_min, h_max
            if dis_cyc:
                self.disruption = cycle(disruption)
            else:
                self.disruption = iter(disruption)
        super(InvertedPendulum, self).__init__(parent=None)

    # Inicjalizacja obrazka
    def init_image(self):
        self.h_scale = self.image_h/(self.h_max-self.h_min)
        self.x_scale = self.image_w/(2*self.x_max)
        self.hor = (self.h_max-10)*self.h_scale
        self.c_w = 16*self.x_scale
        self.c_h = 8*self.h_scale
        self.r = 8
        self.x = self.x0
        self.theta = self.theta0
        self.dx = self.dx0
        self.dtheta = self.dtheta0
        self.setFixedSize(self.image_w, self.image_h)
        self.show()
        self.setWindowTitle("Inverted Pendulum")
        self.update()

    # Rysowanie wahadła i miarki
    def paintEvent(self, e):
        x, x_max, x_scale, theta = self.x, self.x_max, self.x_scale, self.theta
        hor, l, h_scale = self.hor, self.l, self.h_scale
        image_w, c_w, c_h, r, image_h, h_max, h_min = self.image_w, self.c_w, self.c_h, self.r, self.image_h, self.h_max, self.h_min
        painter = QtGui.QPainter(self)
        painter.setPen(pg.mkPen('k', width=2.0*self.h_scale))
        painter.drawLine(0, int(hor), int(image_w), int(hor))
        painter.setPen(pg.mkPen((165, 42, 42), width=2.0*self.x_scale))
        painter.drawLine(int(x_scale*(x+x_max)), int(hor), int(x_scale*(x+x_max-l*sin(theta))), int(hor-h_scale*(l*cos(theta))))
        painter.setPen(pg.mkPen('b'))
        painter.setBrush(pg.mkBrush('b'))
        painter.drawRect(int(x_scale*(x+x_max)-c_w/2), int(hor-c_h/2), int(c_w), int(c_h))
        painter.setPen(pg.mkPen('r'))
        painter.setBrush(pg.mkBrush('r'))
        painter.drawEllipse(int(x_scale*(x+x_max-l*sin(theta)-r/2)), int(hor-h_scale*(l*cos(theta)+r/2)), int(r*x_scale), int(r*h_scale))
        painter.setPen(pg.mkPen('k'))
        for i in np.arange(-x_max, x_max, x_max/10):
            painter.drawText(int((i+x_max)*x_scale), int(image_h-10), str(int(i)))
        for i in np.arange(h_min, h_max, (h_max-h_min)/10):
            painter.drawText(0, int(image_h-(int(i)-h_min)*h_scale), str(int(i)))

    # Rozwiązanie równań mechaniki wahadła
    def solve_equation(self, F):
        l, m, M = self.l, self.m, self.M
        g = 9.81
        a11 = M+m
        a12 = -m*l*cos(self.theta)
        b1 = F-m*l*self.dtheta**2*sin(self.theta)
        a21 = -cos(self.theta)
        a22 = l
        b2 = g*sin(self.theta)
        a = np.array([[a11, a12], [a21, a22]])
        b = np.array([b1, b2])
        sol = np.linalg.solve(a, b)
        return sol[0], sol[1]

    # Scałkowanie numeryczne przyśpieszenia, żeby uzyskać pozostałe parametry układu
    def count_state_params(self, F, dt=0.001):
        ddx, ddtheta = self.solve_equation(F)
        self.dx += ddx*dt
        self.x += self.dx*dt
        self.dtheta += ddtheta*dt
        self.theta += self.dtheta*dt
        self.theta = arctan2(sin(self.theta), cos(self.theta))

    # Uruchomienie symulacji
    # Zmienna sandbox mówi o tym, czy symulacja ma zostać przerwana w przypadku nieudanego sterowania -
    # - to znaczy takiego, które pozwoliło na zbyt duże wychylenia iksa lub na zbyt poziomo położenie wahadła
    def run(self, sandbox, frameskip=20):
        self.sandbox = sandbox
        self.frameskip = frameskip
        self.init_image()
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.single_loop_run)
        timer.start(1)

    # n - krotne obliczenie następnego stanu układu
    # Gdzie n - to frameskip
    def single_loop_run(self):
        for i in range(self.frameskip+1):
            dis=next(self.disruption, 0)

            # print(f"Polozenie wozka : {self.x}, predkosc wozka : {self.dx}, polozenie kątowe wahadla : {self.theta}, predkosc kątowa wahadla : {self.dtheta}")
            # norm_x = self.normalize_inputs(self.x, MIN_X, MAX_X)
            # norm_dx = self.normalize_inputs(self.dx, self.x_min, self.x_max)
            # norm_theta = self.normalize_inputs(self.theta, MIN_THETA, MAX_THETA)
            # norm_dtheta = self.normalize_inputs(self.dtheta, -np.pi, np.pi)

            # control = self.fuzzy_control(norm_x, norm_theta, norm_dx, norm_dtheta)
            control = self.fuzzy_control(self.x, self.theta, self.dx, self.dtheta)
            # print(f"Znormalizowane wartosci: x: {norm_x}, dx: {norm_dx}, theta: {norm_theta}, dtheta: {norm_dtheta}")

            F = dis+control
            self.count_state_params(F)
            if not self.sandbox:
                if self.x < -self.x_max or self.x > self.x_max or np.abs(self.theta) > np.pi/3:
                    exit(1)
        self.update()
        
    def fuzzy_control(self, x, theta, dx, dtheta):
        fuzzy_inputs = self.fuzzify(x, theta, dx, dtheta)
        fuzzy_output = self.inference(fuzzy_inputs)
        output = self.defuzzify(fuzzy_output)
        # return 0
        return output

    def fuzzify(self, x, theta, dx, dtheta):
        fuzzy_values = {
            'x':      self.fuzzify_x(x),
            'theta':  self.fuzzify_theta(theta),
            'dx':     self.fuzzify_dx(dx),
            'dtheta': self.fuzzify_dtheta(dtheta)
        }
        return fuzzy_values

    def fuzzify_x(self, x):
        # print(f"X: {x}")
        fuzzy_x = {
            'negative': self.triangle_membership(x, -self.x_max, -self.x_max/2, 0),
            'zero': self.triangle_membership(x, -self.x_max/2, 0, self.x_max/2),
            'positive': self.triangle_membership(x, 0, self.x_max/2, self.x_max)
        }
        # print(f"Fuzzy x: {fuzzy_x}")
        return fuzzy_x

    def fuzzify_theta(self, theta):
        fuzzy_theta = {
            'negative': self.triangle_membership(theta, -MAX_THETA, -MAX_THETA/2, 0),
            'zero': self.triangle_membership(theta, -MAX_THETA/2, 0, MAX_THETA/2),
            'positive': self.triangle_membership(theta, 0, MAX_THETA/2, MAX_THETA)
        }
        return fuzzy_theta

    def fuzzify_dx(self, dx):
        fuzzy_dx = {
            'negative': self.triangle_membership(dx, MIN_DX, MIN_DX/2, 0),
            'zero': self.triangle_membership(dx, MIN_DX/2, 0, MAX_DX/2),
            'positive': self.triangle_membership(dx, 0, MAX_DX/2, MAX_DX)
        }
        return fuzzy_dx

    def fuzzify_dtheta(self, dtheta):
        # print(f"Dtheta: {dtheta}")
        fuzzy_dtheta = {
            'low': self.triangle_membership(dtheta, MIN_DTHETA, MAX_DTHETA/4, MAX_DTHETA/2),
            'high': self.triangle_membership(dtheta, MAX_DTHETA/2, 3*(MAX_DTHETA/4), MAX_DTHETA)
        }
        return fuzzy_dtheta

    def triangle_membership(self, value, a, b, c):
        if a <= value <= b:
            return (value - a) / (b - a)
        elif b <= value <= c:
            return (c - value) / (c - b)
        else:
            return 0

    def inference(self, fuzzy_inputs):
        rules = [
            (min(fuzzy_inputs['x']['negative'], fuzzy_inputs['dx']['negative']), 'positive_large'),
            (min(fuzzy_inputs['x']['negative'], fuzzy_inputs['dx']['positive']), 'zero'),
            (min(fuzzy_inputs['x']['positive'], fuzzy_inputs['dx']['negative']), 'zero'),
            (min(fuzzy_inputs['x']['positive'], fuzzy_inputs['dx']['positive']), 'negative_large'),

            # (min(fuzzy_inputs['x']['negative'], fuzzy_inputs['theta']['negative']), 'positive_large'),
            # (min(fuzzy_inputs['x']['negative'], fuzzy_inputs['theta']['zero']), 'positive_small'),
            # (min(fuzzy_inputs['x']['negative'], fuzzy_inputs['theta']['positive']), 'positive_small'),
            # (min(fuzzy_inputs['x']['zero'], fuzzy_inputs['theta']['negative']), 'positive_small'),
            # (min(fuzzy_inputs['x']['zero'], fuzzy_inputs['theta']['zero']), 'zero'),
            # (min(fuzzy_inputs['x']['zero'], fuzzy_inputs['theta']['positive']), 'negative_small'),
            # (min(fuzzy_inputs['x']['positive'], fuzzy_inputs['theta']['negative']), 'negative_small'),
            # (min(fuzzy_inputs['x']['positive'], fuzzy_inputs['theta']['zero']), 'negative_small'),
            # (min(fuzzy_inputs['x']['positive'], fuzzy_inputs['theta']['positive']), 'negative_large'),

            (min(fuzzy_inputs['theta']['negative'], fuzzy_inputs['dtheta']['high']), 'positive_small'),
            (min(fuzzy_inputs['theta']['positive'], fuzzy_inputs['dtheta']['low']), 'negative_small'),

            (min(fuzzy_inputs['theta']['negative'], fuzzy_inputs['dtheta']['low']), 'positive_large'),
            (min(fuzzy_inputs['theta']['negative'], fuzzy_inputs['dtheta']['high']), 'positive_small'),
            (min(fuzzy_inputs['theta']['zero'], fuzzy_inputs['dtheta']['low']), 'positive_small'),
            (min(fuzzy_inputs['theta']['zero'], fuzzy_inputs['dtheta']['high']), 'negative_small'),
            (min(fuzzy_inputs['theta']['positive'], fuzzy_inputs['dtheta']['low']), 'negative_small'),
            (min(fuzzy_inputs['theta']['positive'], fuzzy_inputs['dtheta']['high']), 'negative_large')
        ]

        fuzzy_output = {'negative_large': 0, 'negative_small': 0, 'zero': 0, 'positive_small': 0, 'positive_large': 0}
        for rule_strength, output_label in rules:
            fuzzy_output[output_label] = max(fuzzy_output[output_label], rule_strength)

        # print(f"Fuzzy output: {fuzzy_output}")

        return fuzzy_output

    def defuzzify(self, fuzzy_output):
        """ Defuzzification """
        output_values = {'negative_large': -4, 'negative_small': -2, 'zero': 0, 'positive_small': 2, 'positive_large': 4}
        numerator = sum(value * output_values[key] for key, value in fuzzy_output.items())
        denominator = sum(fuzzy_output.values())
        output =  numerator / denominator if denominator != 0 else 0
        print(f"Output: {output}")
        return output

    def normalize_inputs(self, value, min_val, max_val):
        """ Normalization of input values """ # Mistake normalizing 0 to 1
        return (value - min_val) / (max_val - min_val)

if __name__ == '__main__':
    app = QtWidgets.QApplication(argv)
    if len(argv)>1:
        ip = InvertedPendulum(f_name=argv[1])
    else:
        ip = InvertedPendulum(x0=90, dx0=0, theta0=0, dtheta0=0.1, ih=800, iw=1000, h_min=-80, h_max=80)
    ip.run(sandbox=True)
    exit(app.exec_())
