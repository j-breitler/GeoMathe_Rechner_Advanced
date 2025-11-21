#-------------------- GeoMathe_Rechner_Advanced.py -----------#
#-------------------- Author: Ai & Julian Breitler -------------#
#------- Co-Author & Helferlein: Fenja Runfors, Klara Gössler, Christopher Stering---#
#------------------ Viel Spaß beim Verwenden! ----------------#

#-------------------------------------------------------------#
#------------------Benützung auf eigene Gefahr !!! -----------#
#-------------------------------------------------------------#

import math
import tkinter as tk
from tkinter import messagebox, ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# -----------------------------------------------------------
# --- GLOBAL CONSTANTS & PURE CALCULATION FUNCTIONS (HA_combined) ---
# -----------------------------------------------------------

# Define the Gon constant: 400 gon = 2 * pi radians
GON_PER_RADIAN = 200 / math.pi
RADIAN_PER_GON = math.pi / 200

def normalize_azimuth(azimuth_gon):
    """
    Normalizes an azimuth value to the standard geodetic range [0, 400) gon 
    and returns the normalized azimuth and a message if normalization occurred.
    """
    original_az = azimuth_gon
    message = ""
    
    normalized_az = azimuth_gon % 400
    
    # Check if normalization actually occurred for the message
    if normalized_az != original_az:
        # Simplified message generation based on original logic
        if original_az < 0:
            message = f"Modulo operation: {original_az:.4f} + {math.ceil(-original_az / 400) * 400} gon"
        elif original_az >= 400:
            message = f"Modulo operation: {original_az:.4f} - {math.floor(original_az / 400) * 400} gon"
        
    if normalized_az == 400:
        normalized_az = 0.0
        
    return normalized_az, message

def calculate_geodetic_inverse(x1, y1, x2, y2):
    """
    Calculates the distance (S12) and forward/backward azimuths (Az12, Az21).
    """
    dx = x2 - x1
    dy = y2 - y1
    
    if dx == 0 and dy == 0:
        return 0.0, 0.0, 0.0, "", ""
        
    S12 = math.sqrt(dx**2 + dy**2)
    
    Az12_rad = math.atan2(dy, dx)
    Az12_gon_unnormalized = Az12_rad * GON_PER_RADIAN
    
    Az12_gon, Az12_msg = normalize_azimuth(Az12_gon_unnormalized)
        
    Az21_unnormalized = Az12_gon + 200
    Az21_gon, Az21_msg = normalize_azimuth(Az21_unnormalized)

    return S12, Az12_gon, Az21_gon, Az12_msg, Az21_msg

def calculate_geodetic_forward(x1, y1, S12, Az12_gon):
    """
    Calculates the coordinates of P2 (x2, y2) given P1 (x1, y1), 
    distance S12, and forward azimuth Az12.
    """
    Az12_gon_norm, Az12_msg = normalize_azimuth(Az12_gon)
    Az12_rad = Az12_gon_norm * RADIAN_PER_GON
    
    dx = S12 * math.cos(Az12_rad)
    dy = S12 * math.sin(Az12_rad)
    
    x2 = x1 + dx
    y2 = y1 + dy
    
    return x2, y2, Az12_msg

# -----------------------------------------------------------
# --- PURE CALCULATION FUNCTIONS (HWS_berechnung_alle_Winkel) ---
# -----------------------------------------------------------

def calculate_all_angles(a: float, b: float, c: float) -> dict:
    """
    Calculates all three angles (alpha, beta, gamma) of a triangle 
    using the Half-Angle Theorem (Halbwinkelsatz), returning the results in Gon.
    """
    
    # --- Input Validation (Triangle Inequality) ---
    if a + b <= c or a + c <= b or b + c <= a:
        raise ValueError("Ungültiges Dreieck (Dreiecksungleichung verletzt)")

    # --- Calculation of Semi-perimeter and differences ---
    s = (a + b + c) / 2
    s_minus_a = s - a
    s_minus_b = s - b
    s_minus_c = s - c

    # --- General Calculation Helper ---
    def calculate_angle(s_opp_num, s_adj1_num, s_adj2_num):
        """Calculates an angle based on the side differences."""
        zaehler = s_adj1_num * s_adj2_num
        nenner = s * s_opp_num
        
        # Check for non-physical/edge cases before sqrt/atan
        if nenner == 0 or zaehler < 0 or nenner < 0:
             return math.nan, math.nan 

        tan_half = math.sqrt(zaehler / nenner)
        angle_rad = 2 * math.atan(tan_half)
        angle_gon = angle_rad * GON_PER_RADIAN
        return angle_rad, angle_gon
    
    # --- Calculate Angles ---
    alpha_rad, alpha_gon = calculate_angle(s_minus_a, s_minus_b, s_minus_c)
    beta_rad, beta_gon = calculate_angle(s_minus_b, s_minus_a, s_minus_c)
    gamma_rad, gamma_gon = calculate_angle(s_minus_c, s_minus_a, s_minus_b)

    # --- Return Results ---
    return {
        "Semiperimeter (s)": s,
        "s - a": s_minus_a,
        "s - b": s_minus_b,
        "s - c": s_minus_c,
        
        "--- Ergebnisse im Gon-System ---": "", # Separator
        
        "Alpha (α) [Gon]": alpha_gon,
        "Beta (β) [Gon]": beta_gon,
        "Gamma (γ) [Gon]": gamma_gon,
        
        "--- Überprüfung ---": "", # Separator
        
        "Summe (α+β+γ) [Gon]": alpha_gon + beta_gon + gamma_gon
    }

# -----------------------------------------------------------
# --- PURE CALCULATION FUNCTIONS (Numerisch_stabiler_Algorithmus) ---
# -----------------------------------------------------------

def compute_N_from_LMR_gon(xL, yL, xM, yM, xR, yR, alpha_gon, beta_gon):
    """
    Calculates the coordinates of point N using the numerically stable 
    resection algorithm (Pothenot-Snellius).
    """
    
    # 1) Distances s_ML, s_MR from coordinates
    dx_ML = xL - xM
    dy_ML = yL - yM
    dx_MR = xR - xM
    dy_MR = yR - yM

    s_ML = math.hypot(dx_ML, dy_ML)
    s_MR = math.hypot(dx_MR, dy_MR)

    # 2) Convert angles to radians
    alpha = alpha_gon * RADIAN_PER_GON
    beta = beta_gon * RADIAN_PER_GON

    # 3) Azimuths ν_ML, ν_MR (in radians)
    nu_ML = math.atan2(dy_ML, dx_ML)
    nu_MR = math.atan2(dy_MR, dx_MR)

    # 4) a, b
    a = math.sin(alpha) / s_ML
    b = -math.sin(beta) / s_MR

    # 5) Denominator
    den = math.sin(nu_ML - nu_MR + alpha + beta)
    if abs(den) < 1e-12:
        raise ValueError("Denominator ~ 0 (Geometrieproblem / Singularität).")

    # 6) λ and µ
    lam = (a * math.cos(nu_MR - beta) - b * math.cos(nu_ML + alpha)) / den
    mu = (a * math.sin(nu_MR - beta) - b * math.sin(nu_ML + alpha)) / den

    # 7) s_mn^2 (Square of the distance M to N)
    s_mn_sq = 1.0 / (lam**2 + mu**2)

    # 8) Δx_mn, Δy_mn (Coordinate differences M to N)
    dx_mn = lam * s_mn_sq
    dy_mn = mu * s_mn_sq

    # 9) Point N
    x_n = xM + dx_mn
    y_n = yM + dy_mn

    return {
        "s_ML": s_ML,
        "s_MR": s_MR,
        "alpha_gon": alpha_gon,
        "beta_gon": beta_gon,
        "nu_ML_gon": nu_ML * GON_PER_RADIAN,
        "nu_MR_gon": nu_MR * GON_PER_RADIAN,
        "a": a,
        "b": b,
        "Nenner von Formel": den,
        "lambda": lam,
        "mu": mu,
        "s_mn^2": s_mn_sq,
        "dx_mn": dx_mn,
        "dy_mn": dy_mn,
        "x_n": x_n,
        "y_n": y_n,
    }

# -----------------------------------------------------------
# --- MAIN APPLICATION CLASS ---
# -----------------------------------------------------------

class GeodesyApp:
    def __init__(self, root):
        self.root = root
        self.root.title("GeoMathe Rechner - Geodaesie")
        
        # --- MODERNIZATION & FONT UPDATE ---
        # Set a modern theme (e.g., 'clam' or 'alt')
        style = ttk.Style()
        style.theme_use('clam')
        
        # Increase the base font size for better readability
        # Adjust 'TkDefaultFont' and other relevant fonts
        base_font = ('Arial', 14)
        style.configure('.', font=base_font)
        style.configure('TButton', font=('Arial', 14, 'bold'))
        style.configure('TLabel', font=('Arial', 14))
        style.configure('TEntry', font=('Arial', 14))
        style.configure('TCombobox', font=('Arial', 14))
        style.configure('Treeview.Heading', font=('Arial', 14, 'bold'))
        style.configure('Treeview', font=('Arial', 14))
        # -----------------------------------
        
        # Start with an EMPTY coordinate dictionary ---
        self.COORDINATE_DATA = {}
        self.POINT_IDS = sorted(self.COORDINATE_DATA.keys())
        # ---------------------------------------------------------
        
        # Dictionary to hold entry/combobox references for easy access
        self.widgets = {}

        # Create a tabbed interface (Notebook)
        self.notebook = ttk.Notebook(root)
        self.notebook.pack(pady=10, padx=10, expand=True, fill="both")

        # Setup all tabs
        self._setup_point_manager()
        self._setup_inverse_tab()
        self._setup_forward_tab()
        self._setup_half_angle_tab()
        self._setup_n_point_tab()
        
        # Initial GUI synchronization
        self._refresh_comboboxes()
        self.notebook.bind("<<NotebookTabChanged>>", self._on_tab_change)

    # --- Utility Methods ---

    def _on_tab_change(self, event):
        """Called when the user switches tabs to refresh related data."""
        current_tab = self.notebook.tab(self.notebook.select(), "text")
        
        # Refresh coordinate labels if switching to a calculation tab
        if current_tab == '1. Geodaetische Hauptaufgabe':
            self._update_coords_fwd(None)
        elif current_tab == '2. Geodaetische Hauptaufgabe':
            self._update_coords_inv(None)

    def _refresh_comboboxes(self):
        """Updates the list of available points in all relevant Comboboxes."""
        self.POINT_IDS = sorted(self.COORDINATE_DATA.keys())
        
        for key in ['combo_p1_inv', 'combo_p2_inv', 'combo_p1_fwd', 'combo_pL', 'combo_pM', 'combo_pR']:
            if key in self.widgets:
                combo = self.widgets[key]
                combo['values'] = self.POINT_IDS
                # Maintain selection if possible, otherwise set to empty string
                if combo.get() not in self.POINT_IDS:
                    if self.POINT_IDS:
                        combo.set(self.POINT_IDS[0]) # Default to first point if list is not empty
                    else:
                        combo.set("") # Set to empty if list is empty
  

    def _plot_points(self):
        """Creates a scatter plot of all points in COORDINATE_DATA."""
        # Clear the previous plot
        self.ax.clear()
        
        if not self.COORDINATE_DATA:
            self.ax.set_title("No Points to Plot")
            self.ax.set_xlabel("Y [m] (Easting)")
            self.ax.set_ylabel("X [m] (Northing)")
            self.ax.axis('equal')
            self.canvas.draw()
            return

        # Prepare data for plotting
        p_ids = []
        x_coords = []
        y_coords = []
        for p_id, (x, y) in self.COORDINATE_DATA.items():
            p_ids.append(p_id)
            x_coords.append(x)
            y_coords.append(y)
            
        # Create scatter plot - SWAPPED: plot Y on x-axis, X on y-axis
        self.ax.scatter(y_coords, x_coords, label='Points', marker='o', color='tab:blue')
        
        # Add labels to the points - SWAPPED coordinates
        for p_id, x, y in zip(p_ids, x_coords, y_coords):
            # Annotate slightly offset from the point
            self.ax.annotate(p_id, (y, x), textcoords="offset points", xytext=(5,5), ha='left')

        # Configure plot appearance - SWAPPED axis labels
        self.ax.set_title("Point Overview")
        self.ax.set_xlabel("Y [m] (Easting)")
        self.ax.set_ylabel("X [m] (Northing)")
        self.ax.grid(True, linestyle=':', alpha=0.6)
        
        # Ensure equal scaling (CRUCIAL for geodetic visualization)
        self.ax.axis('equal')
        
        # Redraw the canvas
        self.canvas.draw()

    # -----------------------------------------------------------
    # --- Tab 1: Point Management Setup and Logic ---
    # -----------------------------------------------------------

    def _setup_point_manager(self):
        tab_point_manager = ttk.Frame(self.notebook)
        self.notebook.add(tab_point_manager, text='1. Punktverwaltung')

        # --- Configure Grid Layout for Two Halves ---
        tab_point_manager.grid_columnconfigure(0, weight=1) # Left half for input/list
        tab_point_manager.grid_columnconfigure(1, weight=1) # Right half for plot
        tab_point_manager.grid_rowconfigure(0, weight=0) # Input row: fixed size
        tab_point_manager.grid_rowconfigure(1, weight=1) # Treeview/Plot row: expands

        # --- LEFT HALF: INPUT (Fixed widget keys) ---
        
        input_frame = ttk.LabelFrame(tab_point_manager, text="Point Input", padding=(10, 10))
        input_frame.grid(row=0, column=0, sticky=(tk.N, tk.W, tk.E), padx=10, pady=10)

        # Point ID
        ttk.Label(input_frame, text="Point ID:").grid(row=0, column=0, sticky=tk.W, pady=5, padx=5)
        self.widgets['entry_p_id'] = ttk.Entry(input_frame, width=15)
        self.widgets['entry_p_id'].grid(row=0, column=1, sticky=(tk.W, tk.E), pady=5, padx=5)
        
        # X Coordinate
        ttk.Label(input_frame, text="X-Coordinate:").grid(row=1, column=0, sticky=tk.W, pady=5, padx=5)
        self.widgets['entry_x_val'] = ttk.Entry(input_frame, width=15)
        self.widgets['entry_x_val'].grid(row=1, column=1, sticky=(tk.W, tk.E), pady=5, padx=5)
        self.widgets['entry_x_val'].insert(0, "0")
        
        # Y Coordinate
        ttk.Label(input_frame, text="Y-Coordinate:").grid(row=1, column=2, sticky=tk.W, pady=5, padx=5)
        self.widgets['entry_y_val'] = ttk.Entry(input_frame, width=15)
        self.widgets['entry_y_val'].grid(row=1, column=3, sticky=(tk.W, tk.E), pady=5, padx=5)
        self.widgets['entry_y_val'].insert(0, "0")
        
        button_frame = ttk.Frame(input_frame)
        button_frame.grid(row=2, column=0, columnspan=4, pady=10)

        ttk.Button(button_frame, text="Add/Update Point", command=self._add_update_point).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Delete Selected Point", command=self._delete_selected_point).pack(side=tk.LEFT, padx=5)

        # --- LEFT HALF: POINT LIST (Fixed widget keys) ---

        current_points_frame = ttk.LabelFrame(tab_point_manager, text="Current Points", padding=(10, 10))
        current_points_frame.grid(row=1, column=0, sticky=(tk.N, tk.S, tk.W, tk.E), padx=10, pady=10)
        current_points_frame.grid_rowconfigure(0, weight=1)
        current_points_frame.grid_columnconfigure(0, weight=1)
        
        self.widgets['tree_points'] = ttk.Treeview(current_points_frame, columns=('Point ID', 'X', 'Y'), show='tree headings', selectmode='browse')
        self.widgets['tree_points'].heading('#0', text='')
        self.widgets['tree_points'].heading('Point ID', text='Point ID')
        self.widgets['tree_points'].heading('X', text='X-Coordinate')
        self.widgets['tree_points'].heading('Y', text='Y-Coordinate')
        
        self.widgets['tree_points'].column('#0', width=0, stretch=False)
        self.widgets['tree_points'].column('Point ID', width=80, anchor=tk.CENTER)
        self.widgets['tree_points'].column('X', width=120, anchor=tk.CENTER)
        self.widgets['tree_points'].column('Y', width=120, anchor=tk.CENTER)
        
        self.widgets['tree_points'].grid(row=0, column=0, sticky=(tk.N, tk.S, tk.W, tk.E))
        self.widgets['tree_points'].bind('<<TreeviewSelect>>', self._select_point_from_tree)
        
        scrollbar = ttk.Scrollbar(current_points_frame, orient=tk.VERTICAL, command=self.widgets['tree_points'].yview)
        self.widgets['tree_points'].configure(yscrollcommand=scrollbar.set)
        scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))

        # --- RIGHT HALF: PLOT (NEW) ---
        
        plot_frame = ttk.LabelFrame(tab_point_manager, text="Point Overview Plot", padding=(5, 5))
        # Spans both rows in column 1
        plot_frame.grid(row=0, column=1, rowspan=2, sticky=(tk.N, tk.S, tk.W, tk.E), padx=10, pady=10)

        # Initialize Matplotlib Figure and Canvas
        # Use a small figsize, let the frame resize it
        self.fig, self.ax = plt.subplots(figsize=(1, 1), dpi=100) 
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill=tk.BOTH, expand=True)

        # Draw initial plot
        self._plot_points()

    def _add_update_point(self):
        """Reads input fields and adds/updates a point in COORDINATE_DATA."""
        try:
            p_id = self.widgets['entry_p_id'].get().strip()
            x_val = float(self.widgets['entry_x_val'].get())
            y_val = float(self.widgets['entry_y_val'].get())
            
            if not p_id:
                messagebox.showerror("Input Error", "Point ID cannot be empty.")
                return

            self.COORDINATE_DATA[p_id] = (x_val, y_val)
            
            self._load_points_to_treeview()
            self._refresh_comboboxes()
            self._plot_points()
            
        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numbers for X and Y coordinates.")
        except Exception as e:
            messagebox.showerror("Error", f"An unexpected error occurred: {e}")

    def _load_points_to_treeview(self):
        """Clears and reloads all points from COORDINATE_DATA into the Treeview."""
        tree_points = self.widgets['tree_points']
        for i in tree_points.get_children():
            tree_points.delete(i)
            
        for name, (x, y) in self.COORDINATE_DATA.items():
            tree_points.insert("", "end", iid=name, values=(name, f"{x:.4f}", f"{y:.4f}"))

    def _select_point_from_tree(self, event):
        """Loads the selected point from the Treeview back into the input fields for editing."""
        tree_points = self.widgets['tree_points']
        selected_item = tree_points.focus()
        if selected_item:
            values = tree_points.item(selected_item, 'values')
            if values:
                self.widgets['entry_p_id'].delete(0, tk.END)
                self.widgets['entry_x_val'].delete(0, tk.END)
                self.widgets['entry_y_val'].delete(0, tk.END)
                
                self.widgets['entry_p_id'].insert(0, values[0])
                self.widgets['entry_x_val'].insert(0, str(values[1]))
                self.widgets['entry_y_val'].insert(0, str(values[2]))
            
    def _delete_selected_point(self):
        """Deletes the selected point from the Treeview and the data structure."""
        tree_points = self.widgets['tree_points']
        selected_item = tree_points.focus()
        if selected_item:
            p_id = tree_points.item(selected_item, 'values')[0]
            
            if messagebox.askyesno("Confirm Delete", f"Are you sure you want to delete point '{p_id}'?"):
                del self.COORDINATE_DATA[p_id]
                
                self._load_points_to_treeview()
                self._refresh_comboboxes()
                self._plot_points()
                
                self.widgets['entry_p_id'].delete(0, tk.END)
                self.widgets['entry_x_val'].delete(0, tk.END)
                self.widgets['entry_y_val'].delete(0, tk.END)
        else:
            messagebox.showwarning("Selection Error", "Please select a point to delete.")

    # -----------------------------------------------------------
    # --- Tab 2: Inverse Geodetic Problem (2.HA) Setup and Logic ---
    # -----------------------------------------------------------

    def _setup_inverse_tab(self):
        tab_inverse = ttk.Frame(self.notebook)
        self.notebook.add(tab_inverse, text='2. Geodaetische Hauptaufgabe') 

        input_frame_inv = ttk.LabelFrame(tab_inverse, text="Input Point IDs (Pₙ | Pₘ)", padding=(10, 10))
        input_frame_inv.pack(padx=10, pady=10, fill="x")

        # Point 1 (P1 ID)
        ttk.Label(input_frame_inv, text="Point Pₙ:").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.widgets['combo_p1_inv'] = ttk.Combobox(input_frame_inv, values=self.POINT_IDS, state="readonly", width=12)
        self.widgets['combo_p1_inv'].grid(row=0, column=1, padx=5, pady=5)
        self.widgets['combo_p1_inv'].set("") 
        
        # Use ttk.Label for theme support and custom style for coordinates
        self.widgets['coord_x1_inv'] = ttk.Label(input_frame_inv, text=f"X₁: --.--", anchor="w", foreground="#607D8B")
        self.widgets['coord_x1_inv'].grid(row=1, column=0, padx=5, pady=0, sticky="w")
        self.widgets['coord_y1_inv'] = ttk.Label(input_frame_inv, text=f"Y₁: --.--", anchor="w", foreground="#607D8B")
        self.widgets['coord_y1_inv'].grid(row=1, column=1, padx=5, pady=0, sticky="w")

        # Point 2 (P2 ID)
        ttk.Label(input_frame_inv, text="Point Pₘ:").grid(row=0, column=2, padx=5, pady=5, sticky="w")
        self.widgets['combo_p2_inv'] = ttk.Combobox(input_frame_inv, values=self.POINT_IDS, state="readonly", width=12)
        self.widgets['combo_p2_inv'].grid(row=0, column=3, padx=5, pady=5)
        self.widgets['combo_p2_inv'].set("")
        
        # Use ttk.Label for theme support and custom style for coordinates
        self.widgets['coord_x2_inv'] = ttk.Label(input_frame_inv, text=f"X₂: --.--", anchor="w", foreground="#607D8B")
        self.widgets['coord_x2_inv'].grid(row=1, column=2, padx=5, pady=0, sticky="w")
        self.widgets['coord_y2_inv'] = ttk.Label(input_frame_inv, text=f"Y₂: --.--", anchor="w", foreground="#607D8B")
        self.widgets['coord_y2_inv'].grid(row=1, column=3, padx=5, pady=0, sticky="w")

        self.widgets['combo_p1_inv'].bind("<<ComboboxSelected>>", self._update_coords_inv)
        self.widgets['combo_p2_inv'].bind("<<ComboboxSelected>>", self._update_coords_inv)

        # Calculation Button
        calc_button_inv = ttk.Button(tab_inverse, text="Berechne 2.HA (S, νₙₘ, νₘₙ)", command=self._run_inverse_calculation)
        calc_button_inv.pack(padx=10, pady=10, fill="x")

        # Result Frame
        result_frame_inv = ttk.LabelFrame(tab_inverse, text="Ergebnisse (S, νₙₘ, νₘₙ)", padding=(10, 10))
        result_frame_inv.pack(padx=10, pady=10, fill="x")
        
        # Define a style for result labels for emphasis
        style = ttk.Style()
        style.configure('Result.TLabel', foreground="#00796B", font=('Arial', 11, 'bold'))

        ttk.Label(result_frame_inv, text="Strecke Sₙₘ:").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.widgets['result_S12'] = ttk.Label(result_frame_inv, text="--.--", width=15, anchor="w", style='Result.TLabel')
        self.widgets['result_S12'].grid(row=0, column=1, padx=5, pady=5, sticky="w")

        ttk.Label(result_frame_inv, text="Orientierte Richtung νₙₘ:").grid(row=1, column=0, padx=5, pady=5, sticky="w")
        self.widgets['result_Az12'] = ttk.Label(result_frame_inv, text="--.-- gon", width=15, anchor="w", style='Result.TLabel')
        self.widgets['result_Az12'].grid(row=1, column=1, padx=5, pady=5, sticky="w")
        
        self.widgets['msg_Az12'] = tk.Label(result_frame_inv, text="", fg="#D32F2F")
        self.widgets['msg_Az12'].grid(row=1, column=2, padx=5, pady=5, sticky="w")

        ttk.Label(result_frame_inv, text="Orientierte Richtung νₘₙ:").grid(row=2, column=0, padx=5, pady=5, sticky="w")
        self.widgets['result_Az21'] = ttk.Label(result_frame_inv, text="--.-- gon", width=15, anchor="w", style='Result.TLabel')
        self.widgets['result_Az21'].grid(row=2, column=1, padx=5, pady=5, sticky="w")
        
        self.widgets['msg_Az21'] = tk.Label(result_frame_inv, text="", fg="#D32F2F")
        self.widgets['msg_Az21'].grid(row=2, column=2, padx=5, pady=5, sticky="w")

    def _update_coords_inv(self, event):
        """Updates the displayed coordinates based on selected P1 and P2 IDs."""
        p1_id = self.widgets['combo_p1_inv'].get()
        p2_id = self.widgets['combo_p2_inv'].get()
        
        # Define the coordinate display style for consistency
        coord_style = {'foreground': "#607D8B"}

        if p1_id in self.COORDINATE_DATA:
            x1, y1 = self.COORDINATE_DATA[p1_id]
            self.widgets['coord_x1_inv'].config(text=f"X₁: {x1:.4f}", **coord_style)
            self.widgets['coord_y1_inv'].config(text=f"Y₁: {y1:.4f}", **coord_style)
        else:
            self.widgets['coord_x1_inv'].config(text="X₁: --.--", **coord_style)
            self.widgets['coord_y1_inv'].config(text="Y₁: --.--", **coord_style)

        if p2_id in self.COORDINATE_DATA:
            x2, y2 = self.COORDINATE_DATA[p2_id]
            self.widgets['coord_x2_inv'].config(text=f"X₂: {x2:.4f}", **coord_style)
            self.widgets['coord_y2_inv'].config(text=f"Y₂: {y2:.4f}", **coord_style)
        else:
            self.widgets['coord_x2_inv'].config(text="X₂: --.--", **coord_style)
            self.widgets['coord_y2_inv'].config(text="Y₂: --.--", **coord_style)

    def _run_inverse_calculation(self):
        """Runs the Inverse Problem calculation."""
        try:
            p1_id = self.widgets['combo_p1_inv'].get()
            p2_id = self.widgets['combo_p2_inv'].get()
            
            if p1_id not in self.COORDINATE_DATA or p2_id not in self.COORDINATE_DATA:
                messagebox.showerror("Input Error", "Please select valid point IDs that exist in Point Management.")
                return
                
            x1, y1 = self.COORDINATE_DATA[p1_id]
            x2, y2 = self.COORDINATE_DATA[p2_id]
            
            S12, Az12, Az21, Az12_msg, Az21_msg = calculate_geodetic_inverse(x1, y1, x2, y2)
            
            self.widgets['result_S12'].config(text=f"{S12:.4f}")
            self.widgets['result_Az12'].config(text=f"{Az12:.4f} gon")
            self.widgets['result_Az21'].config(text=f"{Az21:.4f} gon")

            self.widgets['msg_Az12'].config(text=f"⤷ {Az12_msg}" if Az12_msg else "")
            self.widgets['msg_Az21'].config(text=f"⤷ {Az21_msg}" if Az21_msg else "")
            
        except Exception as e:
            messagebox.showerror("Calculation Error", f"An unexpected error occurred: {e}")

    # -----------------------------------------------------------
    # --- Tab 3: Forward Geodetic Problem (1.HA) Setup and Logic ---
    # -----------------------------------------------------------
    
    def _setup_forward_tab(self):
        tab_forward = ttk.Frame(self.notebook)
        self.notebook.add(tab_forward, text='1. Geodaetische Hauptaufgabe')

        input_frame_fwd = ttk.LabelFrame(tab_forward, text="Input (P ID | Sₙₘ | νₙₘ)", padding=(10, 10))
        input_frame_fwd.pack(padx=10, pady=10, fill="x")

        # Point 1 (P1 ID)
        ttk.Label(input_frame_fwd, text="Point ID Pₙ:").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.widgets['combo_p1_fwd'] = ttk.Combobox(input_frame_fwd, values=self.POINT_IDS, state="readonly", width=12)
        self.widgets['combo_p1_fwd'].grid(row=0, column=1, padx=5, pady=5)
        self.widgets['combo_p1_fwd'].set("") 
        
        # Use ttk.Label for theme support and custom style for coordinates
        self.widgets['coord_x1_fwd'] = ttk.Label(input_frame_fwd, text=f"X₁: --.--", anchor="w", foreground="#607D8B")
        self.widgets['coord_x1_fwd'].grid(row=1, column=0, padx=5, pady=0, sticky="w")
        self.widgets['coord_y1_fwd'] = ttk.Label(input_frame_fwd, text=f"Y₁: --.--", anchor="w", foreground="#607D8B")
        self.widgets['coord_y1_fwd'].grid(row=1, column=1, padx=5, pady=0, sticky="w")

        self.widgets['combo_p1_fwd'].bind("<<ComboboxSelected>>", self._update_coords_fwd)

        # Distance (S12) and Azimuth (Az12)
        ttk.Label(input_frame_fwd, text="Strecke Sₙₘ:").grid(row=0, column=2, padx=5, pady=5, sticky="w")
        self.widgets['entry_S12_fwd'] = ttk.Entry(input_frame_fwd, width=15)
        self.widgets['entry_S12_fwd'].grid(row=0, column=3, padx=5, pady=5)
        self.widgets['entry_S12_fwd'].insert(0, "0") 

        ttk.Label(input_frame_fwd, text="Orientierte Richtung νₙₘ (gon):").grid(row=1, column=2, padx=5, pady=5, sticky="w")
        self.widgets['entry_Az12_fwd'] = ttk.Entry(input_frame_fwd, width=15)
        self.widgets['entry_Az12_fwd'].grid(row=1, column=3, padx=5, pady=5)
        self.widgets['entry_Az12_fwd'].insert(0, "0") 

        # Calculation Button
        calc_button_fwd = ttk.Button(tab_forward, text="Berechne 1.HA (X_neu, Y_neu)", command=self._run_forward_calculation)
        calc_button_fwd.pack(padx=10, pady=10, fill="x")

        # Result Frame
        result_frame_fwd = ttk.LabelFrame(tab_forward, text="Ergebnisse (P: X_neu, Y_neu)", padding=(10, 10))
        result_frame_fwd.pack(padx=10, pady=10, fill="x")

        ttk.Label(result_frame_fwd, text="X_neu (Northing):").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.widgets['result_x₂'] = ttk.Label(result_frame_fwd, text="--.--", width=15, anchor="w", style='Result.TLabel')
        self.widgets['result_x₂'].grid(row=0, column=1, padx=5, pady=5, sticky="w")

        ttk.Label(result_frame_fwd, text="Y_neu (Easting):").grid(row=1, column=0, padx=5, pady=5, sticky="w")
        self.widgets['result_y₂'] = ttk.Label(result_frame_fwd, text="--.--", width=15, anchor="w", style='Result.TLabel')
        self.widgets['result_y₂'].grid(row=1, column=1, padx=5, pady=5, sticky="w")

        self.widgets['msg_Az12_fwd'] = tk.Label(result_frame_fwd, text="", fg="#D32F2F")
        self.widgets['msg_Az12_fwd'].grid(row=2, column=0, columnspan=2, padx=5, pady=5, sticky="w")

    def _update_coords_fwd(self, event):
        """Updates the displayed P1 coordinates based on selected P1 ID."""
        p1_id = self.widgets['combo_p1_fwd'].get()
        
        # Define the coordinate display style for consistency
        coord_style = {'foreground': "#607D8B"}

        if p1_id in self.COORDINATE_DATA:
            x1, y1 = self.COORDINATE_DATA[p1_id]
            self.widgets['coord_x1_fwd'].config(text=f"X₁: {x1:.4f}", **coord_style)
            self.widgets['coord_y1_fwd'].config(text=f"Y₁: {y1:.4f}", **coord_style)
        else:
            self.widgets['coord_x1_fwd'].config(text="X₁: --.--", **coord_style)
            self.widgets['coord_y1_fwd'].config(text="Y₁: --.--", **coord_style)

    def _run_forward_calculation(self):
        """Runs the Forward Problem calculation."""
        try:
            p1_id = self.widgets['combo_p1_fwd'].get()
            
            S12 = float(self.widgets['entry_S12_fwd'].get())
            Az12_gon = float(self.widgets['entry_Az12_fwd'].get())
            
            if p1_id not in self.COORDINATE_DATA:
                messagebox.showerror("Input Error", "Please select a valid point ID for P1 that exists in Point Management.")
                return

            x1, y1 = self.COORDINATE_DATA[p1_id]
            
            if S12 < 0:
                raise ValueError("Distance (S) must be non-negative.")
                
            x2, y2, Az12_msg = calculate_geodetic_forward(x1, y1, S12, Az12_gon)
            
            self.widgets['result_x₂'].config(text=f"{x2:.4f} m")
            self.widgets['result_y₂'].config(text=f"{y2:.4f} m")

            self.widgets['msg_Az12_fwd'].config(text=f"⤷ Normalized Az: {Az12_msg}" if Az12_msg else "")
            
        except ValueError as ve:
            messagebox.showerror("Input Error", f"Please enter valid numbers/selections: {ve}")
        except Exception as e:
            messagebox.showerror("Calculation Error", f"An unexpected error occurred: {e}")

    # -----------------------------------------------------------
    # --- Tab 4: Half-Angle Theorem (HWS) Setup and Logic ---
    # -----------------------------------------------------------

    def _setup_half_angle_tab(self):
        tab_hws = ttk.Frame(self.notebook)
        self.notebook.add(tab_hws, text='Halbwinkelsatz (HWS)')

        # Use a main PanedWindow to separate Input/Plot from Output
        main_pane = ttk.PanedWindow(tab_hws, orient=tk.VERTICAL)
        main_pane.pack(fill="both", expand=True, padx=10, pady=10)

        # --- TOP FRAME: Input and Plot ---
        top_frame = ttk.Frame(main_pane)
        main_pane.add(top_frame, weight=1) # Allow this frame to expand vertically

        # Left side: Input fields
        input_frame = ttk.LabelFrame(top_frame, text="Eingabe Seitenlängen", padding=10)
        input_frame.pack(side="left", fill="y", padx=(0, 10))

        ttk.Label(input_frame, text="Für Komma bitte Punkt '.' verwenden").grid(
            row=0, column=0, columnspan=2, sticky="w", pady=(0, 10)
        )

        # Input fields
        ttk.Label(input_frame, text="Seite a:").grid(row=1, column=0, sticky="w", padx=5, pady=5)
        self.widgets['entry_a_hws'] = ttk.Entry(input_frame, width=15)
        self.widgets['entry_a_hws'].grid(row=1, column=1, padx=5, pady=5)
        self.widgets['entry_a_hws'].insert(0, "0")

        ttk.Label(input_frame, text="Seite b:").grid(row=2, column=0, sticky="w", padx=5, pady=5)
        self.widgets['entry_b_hws'] = ttk.Entry(input_frame, width=15)
        self.widgets['entry_b_hws'].grid(row=2, column=1, padx=5, pady=5)
        self.widgets['entry_b_hws'].insert(0, "0")

        ttk.Label(input_frame, text="Seite c:").grid(row=3, column=0, sticky="w", padx=5, pady=5)
        self.widgets['entry_c_hws'] = ttk.Entry(input_frame, width=15)
        self.widgets['entry_c_hws'].grid(row=3, column=1, padx=5, pady=5)
        self.widgets['entry_c_hws'].insert(0, "0")

        # Button
        calc_btn = ttk.Button(input_frame, text="Berechnen", command=self._run_half_angle_calculation)
        calc_btn.grid(row=4, column=0, columnspan=2, pady=(10, 0), sticky="ew")

        # Right side: Triangle Plot
        plot_frame = ttk.LabelFrame(top_frame, text="Allgemeines Dreieck", padding=5)
        plot_frame.pack(side="right", fill="both", expand=True)

        # Matplotlib canvas container
        self.widgets['canvas_hws_container'] = ttk.Frame(plot_frame)
        self.widgets['canvas_hws_container'].pack(fill="both", expand=True)
        
        # Initial plot
        self._create_triangle_plot()

        # --- BOTTOM FRAME: Output Text ---
        output_frame_container = ttk.LabelFrame(main_pane, text="Ergebnisse", padding=10)
        main_pane.add(output_frame_container, weight=1) # Allow this frame to expand vertically

        self.widgets['output_hws'] = tk.Text(output_frame_container, height=15, width=50, font=('Courier', 10), wrap="word", borderwidth=2, relief="groove")
        scrollbar = ttk.Scrollbar(output_frame_container, orient="vertical", command=self.widgets['output_hws'].yview)
        self.widgets['output_hws'].config(yscrollcommand=scrollbar.set)
        
        self.widgets['output_hws'].pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Add tagging for better readability (optional, but nice)
        self.widgets['output_hws'].tag_configure('separator', font=('Courier', 10, 'bold'), foreground='#4CAF50')


    def _create_triangle_plot(self):
        """
        Creates and embeds a static matplotlib plot of a generic triangle 
        with standard geodetic labels (A, B, C; a, b, c; alpha, beta, gamma).
        This plot is based on the general representation provided.
        """
        
        # Clear existing canvas if it exists (for updates or initial setup)
        if 'canvas_hws' in self.widgets:
            self.widgets['canvas_hws'].get_tk_widget().destroy()
            plt.close(self.widgets['fig_hws']) # Close the old figure
        
        # --- Define vertices for a generic non-right, non-equilateral triangle ---
        # Vertices (A, B, C)
        P = [
            (0, 0),    # A (Angle alpha, side c, b)
            (10, 0),   # B (Angle beta, side c, a)
            (3, 7)     # C (Angle gamma, side a, b)
        ]
        
        # Unpack coordinates
        x_coords = [p[0] for p in P] + [P[0][0]]
        y_coords = [p[1] for p in P] + [P[0][1]]

        # --- Matplotlib Figure and Axis ---
        fig, ax = plt.subplots(figsize=(4, 3))
        self.widgets['fig_hws'] = fig # Store figure reference to close later
        
        ax.plot(x_coords, y_coords, color='black', linewidth=1.5)
        
        # --- Labels for Vertices (A, B, C) ---
        ax.text(P[0][0] - 0.5, P[0][1] - 0.5, 'A', fontsize=12, fontweight='bold')
        ax.text(P[1][0] + 0.1, P[1][1] - 0.5, 'B', fontsize=12, fontweight='bold')
        ax.text(P[2][0] - 0.2, P[2][1] + 0.3, 'C', fontsize=12, fontweight='bold')
        
        # --- Labels for Sides (a, b, c) ---
        # Side a (opposite A, between B and C)
        ax.text((P[1][0] + P[2][0]) / 2 + 2, (P[1][1] + P[2][1]) / 2 - 0.8, 
                r'$a$', fontsize=12, color='darkred', ha='center', va='center')
        
        # Side b (opposite B, between A and C)
        ax.text((P[0][0] + P[2][0]) / 2 - 0.5, (P[0][1] + P[2][1]) / 2 + 0.1, 
                r'$b$', fontsize=12, color='darkred', ha='center', va='center')

        # Side c (opposite C, between A and B)
        ax.text((P[0][0] + P[1][0]) / 2, (P[0][1] + P[1][1]) / 2 - 0.5, 
                r'$c$', fontsize=12, color='darkred', ha='center', va='center')

        # --- Labels for Angles (alpha, beta, gamma) using LaTeX for Greek letters ---
        # Alpha (at A)
        ax.text(P[0][0] + 0.5, P[0][1] + 0.3, r'$\alpha$', fontsize=12, color='darkblue')

        # Beta (at B)
        ax.text(P[1][0] - 1.0, P[1][1] + 0.3, r'$\beta$', fontsize=12, color='darkblue')

        # Gamma (at C)
        ax.text(P[2][0] + 0.3, P[2][1] - 1.0, r'$\gamma$', fontsize=12, color='darkblue')

        # --- Final Plot Styling ---
        ax.set_aspect('equal', adjustable='box')
        ax.axis('off') # Hide axes, ticks, and borders
        fig.tight_layout(pad=0.5)

        # --- Embed into Tkinter ---
        canvas = FigureCanvasTkAgg(fig, master=self.widgets['canvas_hws_container'])
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack(fill=tk.BOTH, expand=True)
        self.widgets['canvas_hws'] = canvas # Store canvas reference
        canvas.draw()
        
    def _run_half_angle_calculation(self):
        try:
            a = float(self.widgets['entry_a_hws'].get())
            b = float(self.widgets['entry_b_hws'].get())
            c = float(self.widgets['entry_c_hws'].get())
            
            if a <= 0 or b <= 0 or c <= 0:
                raise ValueError("Side lengths must be positive.")

            result = calculate_all_angles(a, b, c)

            self.widgets['output_hws'].delete(1.0, tk.END)
            for key, value in result.items():
                if key.startswith("---"):
                    self.widgets['output_hws'].insert(tk.END, f"\n{key}\n", 'separator')
                else:
                    self.widgets['output_hws'].insert(tk.END, f"{key:30}: {value:.4f}\n")
            
            # Add tagging for better readability (optional, but nice)
            self.widgets['output_hws'].tag_configure('separator', font=('Courier', 10, 'bold'), foreground='#4CAF50')

        except ValueError as e:
            messagebox.showerror("Fehler", str(e))
        except Exception:
            messagebox.showerror("Fehler", "Ungültige Eingabe oder unerwarteter Fehler!")
            

    # -----------------------------------------------------------
    # --- Tab 5: N-Point Calculation Setup and Logic ---
    # -----------------------------------------------------------

    def _update_coords_n_point(self, event):
        """Updates the displayed coordinates for L, M, and R based on selected IDs."""
        points = [
            ('combo_pL', 'coord_xL_npoint', 'coord_yL_npoint', 'L'),
            ('combo_pM', 'coord_xM_npoint', 'coord_yM_npoint', 'M'),
            ('combo_pR', 'coord_xR_npoint', 'coord_yR_npoint', 'R'),
        ]
        coord_style = {'foreground': "#607D8B"}
        
        for combo_key, x_key, y_key, p_char in points:
            # Check if the combobox key exists before trying to get its value
            if combo_key not in self.widgets:
                continue 

            p_id = self.widgets[combo_key].get()
            if p_id in self.COORDINATE_DATA:
                x, y = self.COORDINATE_DATA[p_id]
                self.widgets[x_key].config(text=f"X_{p_char}: {x:.4f}", **coord_style)
                self.widgets[y_key].config(text=f"Y_{p_char}: {y:.4f}", **coord_style)
            else:
                self.widgets[x_key].config(text=f"X_{p_char}: --.--", **coord_style)
                self.widgets[y_key].config(text=f"Y_{p_char}: --.--", **coord_style)

    def _setup_n_point_tab(self):
        tab_n_point = ttk.Frame(self.notebook)
        self.notebook.add(tab_n_point, text='Numerisch-stabiler Algorithmus')

        main_frame = ttk.Frame(tab_n_point, padding="10")
        main_frame.pack(fill="both", expand=True)
        
        row = 0

        # Coordinate inputs frame
        coord_frame = ttk.LabelFrame(main_frame, text="Known Points (L, M, R)", padding=(10, 10))
        coord_frame.grid(row=row, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=10)
        row += 1

        r = 0
        
        # Point L
        ttk.Label(coord_frame, text="Punkt L:").grid(row=r, column=0, sticky=tk.W, pady=5, padx=5)
        self.widgets['combo_pL'] = ttk.Combobox(coord_frame, values=self.POINT_IDS, state="readonly", width=12)
        self.widgets['combo_pL'].grid(row=r, column=1, sticky=(tk.W, tk.E), pady=5, padx=5)
        self.widgets['combo_pL'].bind("<<ComboboxSelected>>", self._update_coords_n_point)
        self.widgets['coord_xL_npoint'] = ttk.Label(coord_frame, text=f"Lx: --.--", anchor="w", foreground="#607D8B")
        self.widgets['coord_xL_npoint'].grid(row=r, column=2, padx=5, pady=0, sticky="w")
        self.widgets['coord_yL_npoint'] = ttk.Label(coord_frame, text=f"Ly: --.--", anchor="w", foreground="#607D8B")
        self.widgets['coord_yL_npoint'].grid(row=r, column=3, padx=5, pady=0, sticky="w")
        r += 1 

        # Point M
        ttk.Label(coord_frame, text="Punkt M:").grid(row=r, column=0, sticky=tk.W, pady=5, padx=5)
        self.widgets['combo_pM'] = ttk.Combobox(coord_frame, values=self.POINT_IDS, state="readonly", width=12)
        self.widgets['combo_pM'].grid(row=r, column=1, sticky=(tk.W, tk.E), pady=5, padx=5)
        self.widgets['combo_pM'].bind("<<ComboboxSelected>>", self._update_coords_n_point)
        self.widgets['coord_xM_npoint'] = ttk.Label(coord_frame, text=f"Mx: --.--", anchor="w", foreground="#607D8B")
        self.widgets['coord_xM_npoint'].grid(row=r, column=2, padx=5, pady=0, sticky="w")
        self.widgets['coord_yM_npoint'] = ttk.Label(coord_frame, text=f"My: --.--", anchor="w", foreground="#607D8B")
        self.widgets['coord_yM_npoint'].grid(row=r, column=3, padx=5, pady=0, sticky="w")
        r += 1 

        # Point R
        ttk.Label(coord_frame, text="Punkt R:").grid(row=r, column=0, sticky=tk.W, pady=5, padx=5)
        self.widgets['combo_pR'] = ttk.Combobox(coord_frame, values=self.POINT_IDS, state="readonly", width=12)
        self.widgets['combo_pR'].grid(row=r, column=1, sticky=(tk.W, tk.E), pady=5, padx=5)
        self.widgets['combo_pR'].bind("<<ComboboxSelected>>", self._update_coords_n_point)
        self.widgets['coord_xR_npoint'] = ttk.Label(coord_frame, text=f"Rx: --.--", anchor="w", foreground="#607D8B")
        self.widgets['coord_xR_npoint'].grid(row=r, column=2, padx=5, pady=0, sticky="w")
        self.widgets['coord_yR_npoint'] = ttk.Label(coord_frame, text=f"Ry: --.--", anchor="w", foreground="#607D8B")
        self.widgets['coord_yR_npoint'].grid(row=r, column=3, padx=5, pady=0, sticky="w")
        r += 1

        # Angle inputs frame 
        angle_frame = ttk.LabelFrame(main_frame, text="Angles measured at N (in Gon)", padding=(10, 10))
        angle_frame.grid(row=row, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=10)
        row += 1

        r = 0
        # Alpha
        ttk.Label(angle_frame, text="α (Gon):").grid(row=r, column=0, sticky=tk.W, pady=5, padx=5)
        self.widgets['alpha_gon'] = ttk.Entry(angle_frame, width=15)
        self.widgets['alpha_gon'].grid(row=r, column=1, sticky=(tk.W, tk.E), pady=5, padx=5)
        self.widgets['alpha_gon'].insert(0, "0")
        
        # Beta (same row)
        ttk.Label(angle_frame, text="β(Gon):").grid(row=r, column=2, sticky=tk.W, pady=5, padx=5)
        self.widgets['beta_gon'] = ttk.Entry(angle_frame, width=15)
        self.widgets['beta_gon'].grid(row=r, column=3, sticky=(tk.W, tk.E), pady=5, padx=5)
        self.widgets['beta_gon'].insert(0, "0")
        
        # Calculate button
        calc_button = ttk.Button(main_frame, text="Berechne Punkt N (Rückwärtsschnitt)", command=self._run_n_point_calculation)
        calc_button.grid(row=row, column=0, columnspan=2, pady=20, sticky="ew")
        row += 1
        
        # Results frame
        results_frame = ttk.LabelFrame(main_frame, text="Ergebnisse", padding=(10, 10))
        results_frame.grid(row=row, column=0, columnspan=2, sticky=(tk.W, tk.E, tk.N, tk.S), pady=10)
        main_frame.grid_rowconfigure(row, weight=1)
        
        # Results text widget with scrollbar (using tk.Text for better formatting control)
        self.widgets['results_text_n_point'] = tk.Text(results_frame, height=15, width=60, font=('Courier', 10), wrap="word", borderwidth=2, relief="groove")
        scrollbar = ttk.Scrollbar(results_frame, orient="vertical", command=self.widgets['results_text_n_point'].yview)
        self.widgets['results_text_n_point'].configure(yscrollcommand=scrollbar.set)
        
        self.widgets['results_text_n_point'].pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")


    def _run_n_point_calculation(self):
        """Get values from inputs and perform calculation for N-Point."""
        try:
            # Get point IDs from comboboxes
            pL_id = self.widgets['combo_pL'].get()
            pM_id = self.widgets['combo_pM'].get()
            pR_id = self.widgets['combo_pR'].get()
            
            # Validate that points exist
            if pL_id not in self.COORDINATE_DATA or pM_id not in self.COORDINATE_DATA or pR_id not in self.COORDINATE_DATA:
                messagebox.showerror("Input Error", "Please select valid point IDs for L, M, and R.")
                return
            
            # Get coordinates from dictionary
            xL, yL = self.COORDINATE_DATA[pL_id]
            xM, yM = self.COORDINATE_DATA[pM_id]
            xR, yR = self.COORDINATE_DATA[pR_id]
            
            # Get angles from entry fields
            alpha_gon = float(self.widgets['alpha_gon'].get())
            beta_gon = float(self.widgets['beta_gon'].get())
            
            # Perform calculation
            result = compute_N_from_LMR_gon(
                xL, yL, xM, yM, xR, yR, alpha_gon, beta_gon
            )
            
            # Display results
            self._display_n_point_results(result)
            
        except ValueError as e:
            if "could not convert" in str(e):
                messagebox.showerror("Fehler", "Bitte geben Sie gültige Zahlenwerte ein!")
            else:
                messagebox.showerror("Fehler", str(e))
        except Exception as e:
            messagebox.showerror("Fehler", f"Ein Fehler ist aufgetreten:\n{str(e)}")
    
    def _display_n_point_results(self, result):
        """Display calculation results in the text widget."""
        self.widgets['results_text_n_point'].delete(1.0, tk.END)
        
        # Add tagging for better readability
        self.widgets['results_text_n_point'].tag_configure('header', font=('Courier', 11, 'bold'), foreground='#004D40')
        self.widgets['results_text_n_point'].tag_configure('result', font=('Courier', 11, 'bold'), foreground='#C2185B')
        self.widgets['results_text_n_point'].tag_configure('separator', font=('Courier', 10, 'bold'), foreground='#4CAF50')
        
        output = "=" * 50 + "\n"
        output += "BERECHNUNGSERGEBNISSE\n"
        output += "=" * 50 + "\n\n"
        self.widgets['results_text_n_point'].insert(tk.END, output, 'header')
        
        self.widgets['results_text_n_point'].insert(tk.END, "Koordinaten von Punkt N:\n")
        
        self.widgets['results_text_n_point'].insert(tk.END, f"  x_n = ", '')
        self.widgets['results_text_n_point'].insert(tk.END, f"{result['x_n']:.4f}\n", 'result')
        
        self.widgets['results_text_n_point'].insert(tk.END, f"  y_n = ", '')
        self.widgets['results_text_n_point'].insert(tk.END, f"{result['y_n']:.4f}\n\n", 'result')
        
        
        output_sep = "-" * 50 + "\n"
        output_sep += "Zwischenwerte:\n"
        output_sep += "-" * 50 + "\n"
        self.widgets['results_text_n_point'].insert(tk.END, output_sep, 'separator')
        
        for key, value in result.items():
            if key not in ['x_n', 'y_n']:
                if isinstance(value, float):
                    self.widgets['results_text_n_point'].insert(tk.END, f"{key:15s} = {value:12.6f}\n")
                else:
                    self.widgets['results_text_n_point'].insert(tk.END, f"{key:15s} = {value}\n")


if __name__ == "__main__":
    root = tk.Tk()
    app = GeodesyApp(root)
    root.mainloop()