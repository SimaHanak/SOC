import numpy as np
import numpy.ma as ma
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# --- Define Parameters class ---
class Params:
    def __init__(self):
        self.M = 1
        self.J = 0.3
        self.M2 = -0.1
        self.S3 = 0.05
        self.M4 = 0.01
        self.alpha = 1
        self.beta = 1
        self.gamma = 1
    
    def kerr_params(self):
        j = self.J / self.M**2
        self.alpha = self.M2 / (self.M**3 * j**2)
        self.beta = self.S3 / (self.M**4 * j**3)
        self.gamma = self.M4 / (self.M**5 * j**4)

p = Params()
p.kerr_params()

# --- Define metric-related functions ---
def A(r, z):
    return (8*r**2*z**2*(24*p.J**2*p.M + 17*p.M**2*p.M2 + 21*p.M4) +
            r**4*(-10*p.J**2*p.M + 7*p.M**5 + 32*p.M2*p.M**2 - 21*p.M4) +
            8*z**4*(20*p.J**2*p.M - 7*p.M**5 - 22*p.M2*p.M**2 - 7*p.M4))

def B(r, z):
    return (r**4*(10*p.J**2*p.M**2 + 10*p.M2*p.M**3 + 21*p.M4*p.M + 7*p.M2**2) +
            4*z**4*(-40*p.J**2*p.M**2 - 14*p.J*p.S3 + 7*p.M**6 + 30*p.M2*p.M**3 + 14*p.M4*p.M + 7*p.M2**2) -
            4*r**2*z**2*(27*p.J**2*p.M**2 - 21*p.J*p.S3 + 7*p.M**6 + 48*p.M2*p.M**3 + 42*p.M4*p.M + 7*p.M2**2))

def H(r, z):
    return (4*r**2*z**2*(p.J*(p.M2 - 2*p.M**3) - 3*p.M*p.S3) +
            r**4*(p.J*p.M2 + 3*p.M*p.S3))

def G(r, z):
    return (r**2*(p.J**3*(-(r**4 + 8*z**4 - 12*r**2*z**2)) +
                  p.J*p.M*((p.M**3 + 2*p.M2)*r**4 - 8*(3*p.M**3 + 2*p.M2)*z**4 + 4*(p.M**3 + 10*p.M2)*r**2*z**2) +
                  p.M**2*p.S3*(3*r**4 - 40*z**4 + 12*r**2*z**2)))

def F(r, z):
    return (r**4*(p.S3 - p.J*p.M**2) - 4*r**2*z**2*(p.J*p.M**2 + p.S3))

def f(r, z):
    return (1 - (2*p.M)/np.sqrt(r**2 + z**2) + (2*p.M**2)/(r**2 + z**2) +
            ((p.M2 - p.M**3)*r**2 - 2*(p.M**3 + p.M2)*z**2)/((r**2 + z**2)**(5/2)) +
            (2*z**2*(-p.J**2 + p.M**4 + 2*p.M2*p.M) - 2*p.M*p.M2*r**2)/((r**2 + z**2)**3) +
            A(r, z)/(28*(r**2 + z**2)**(9/2)) + B(r, z)/(14*(r**2 + z**2)**5))

def omega(r, z):
    return (-2*p.J*r**2/(r**2 + z**2)**(3/2) -
            2*p.J*p.M*r**2/(r**2 + z**2)**2 +
            F(r, z)/(r**2 + z**2)**(7/2) +
            H(r, z)/(2*(r**2 + z**2)**4) +
            G(r, z)/(4*(r**2 + z**2)**(11/2)))

def gamma(r, z):
    return ((r**2*(p.J**2*(r**2 - 8*z**2) + p.M*(p.M**3 + 3*p.M2)*(r**2 - 4*z**2)))/(4*(r**2 + z**2)**4) - 
            p.M**2*r**2/(2*(r**2 + z**2)**2))

def g_tt(r, z): return -f(r, z)
def g_tf(r, z): return omega(r, z)*f(r, z)
def g_ff(r, z): return -f(r, z)*omega(r, z)**2 + r**2/f(r, z)
def g_rr(r, z): return np.exp(2*gamma(r, z))/f(r, z)
def g_zz(r, z): return np.exp(2*gamma(r, z))/f(r, z)

def Det(r, z): return g_tt(r, z)*g_ff(r, z) - g_tf(r, z)**2

def V_eff(r, z, E, L_z):
    return 1/g_rr(r, z)*(1 + (g_ff(r, z)*E**2 + g_tt(r, z)*L_z**2 + 2*g_tf(r, z)*E*L_z)/(Det(r, z)))

# --- Grids ---
r = np.linspace(1, 15, 150)
z = np.linspace(-8, 8, 150)
r_grid, z_grid = np.meshgrid(r, z)

# --- Start Dash app ---
app = dash.Dash(__name__)

# --- Layout ---
app.layout = html.Div([
    html.H1("Effective Potential Explorer", style={'text-align': 'center'}),

    html.Div([
        html.Label('Energy E'),
        dcc.Slider(id='E', min=0, max=1, step=0.005, value=0.5, marks={0: '0', 1: '1'}, tooltip={'always_visible': True}),

        html.Label('Angular Momentum L_z'),
        dcc.Slider(id='Lz', min=0, max=5, step=0.01, value=3, marks={0: '0', 5: '5'}, tooltip={'always_visible': True}),
        
        html.Label('Spin J'),
        dcc.Slider(id='J', min=0, max=1, step=0.01, value=p.J, marks={0: '0', 1: '1'}, tooltip={'always_visible': True}),

        html.Label('Mass Quadrupole M2'),
        dcc.Slider(id='M2', min=-1, max=1, step=0.01, value=p.M2, marks={-1: '-1', 1: '1'}, tooltip={'always_visible': True}),

        html.Label('Spin Octupole S3'),
        dcc.Slider(id='S3', min=-1, max=1, step=0.01, value=p.S3, marks={-1: '-1', 1: '1'}, tooltip={'always_visible': True}),
        
        html.Label('Mass Hexadecapole M4'),
        dcc.Slider(id='M4', min=-1, max=1, step=0.01, value=p.M4, marks={-1: '-1', 1: '1'}, tooltip={'always_visible': True}),
    ], style={'padding': 10, 'flex': 1}),

    html.Div([
        html.P(f"Current Parameters in Kerr-like fashion: alpha={p.alpha}, beta={p.beta}, gamma={p.gamma}")
    ]),

    html.Div([
        dcc.Graph(id='plot-area')
    ], style={'width': '100%', 'display': 'inline-block'}),
])

# --- Callback ---
@app.callback(
    Output('plot-area', 'figure'),
    [Input('E', 'value'),
     Input('Lz', 'value'),
     Input('J', 'value'),
     Input('M2', 'value'),
     Input('S3', 'value'),
     Input('M4', 'value')]
)
def update_graph(E, Lz, J, M2, S3, M4):
    # Update parameters
    p.J = J
    p.M2 = M2
    p.S3 = S3
    p.M4 = M4

    p.kerr_params()

    Veff_line = V_eff(r, 0, E, Lz)
    Veff_surface = V_eff(r_grid, z_grid, E, Lz)

    fig = make_subplots(
            rows=2, cols=2,
            specs=[[{"type": "xy"}, {"type": "xy"}],
                [{"type": "surface"}, None]],
            subplot_titles=("Effective Potential V_eff(r, z=0)", "CZV", "Effective Potential V_eff(r, z)")
        )
    
    # --- 2D Plot ---
    fig.add_trace(go.Scatter(x=r, y=Veff_line, mode='lines', line=dict(color='blue')), row=1, col=1)
    fig.update_layout(title="Effective Potential V_eff(r, z=0)", xaxis_title="r", yaxis_title="V_eff")

    # --- 3D Plot ---
    fig.add_trace(go.Surface(z=Veff_surface, x=r, y=z, colorscale='Viridis', cmin=-0.2, cmax=0.2), row=2, col=1)
    fig.update_layout(title="Effective Potential V_eff(r, z)", scene=dict(
        xaxis_title='r',
        yaxis_title='z',
        zaxis_title='V_eff',
        zaxis_range=[-0.2, 0.2]
    ))

    fig.add_trace(go.Contour(
                z=Veff_surface, x=r, y=z, 
                contours_coloring='lines',  # draw only contour lines (no fill)
                line=dict(color='red', width=2),  # red contour line
                showscale=False,     # hide the color scale bar
                contours=dict(
                    start=0,         # start at V_eff = 0
                    end=0,           # end at V_eff = 0 (same as start for a single level)
                    size=1           # step size (irrelevant when start=end, but required)
                )
            ), row=1, col=2)
    fig.update_layout(title="CZV", scene=dict(
        xaxis_title='r',
        yaxis_title='z',
    ))

    fig.update_layout(height=900, width=1200)
    return fig

# --- Run server ---
if __name__ == '__main__':
    app.run(debug=True)