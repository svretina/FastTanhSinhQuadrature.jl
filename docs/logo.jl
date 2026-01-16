using Luxor

# Definitions from the package
@inline function ordinate(t)
    return tanh(π / 2 * sinh(t))
end
@inline function weight(t)
    tmp = cosh(π / 2 * sinh(t))
    return ((π / 2) * cosh(t)) / (tmp * tmp)
end

# Colors
const C_PURPLE = "#9558B2"
const C_GREEN  = "#389826"
const C_RED    = "#CB3C33"

function draw_logo_centered(filename)
    Drawing(600, 500, filename)
    origin()
    # background("white") # Transparent background
    
    # Scales
    x_scale = 100  
    y_scale = 100  
    
    # Calculate vertical shift
    # We want the S-curve to sit in the vertical middle of the bell curve.
    # The bell curve goes from y=0 to y=weight(0).
    bell_peak = weight(0.0)
    # Luxor y-axis is inverted (up is negative), so we shift "up" by subtracting.
    vertical_offset = (bell_peak / 2) * y_scale
    
    # -----------------------------------------------------
    # Layer 1: The Weight Function (Purple)
    # -----------------------------------------------------
    sethue(C_PURPLE)
    
    t_start, t_end = -3.0, 3.0
    t_range = range(t_start, t_end, length=120)
    
    # Draw Bell (Base remains at y=0 relative to origin)
    move(Point(t_range[1] * x_scale, 0)) 
    for t in t_range
        y = weight(t)
        line(Point(t * x_scale, -y * y_scale))
    end
    line(Point(t_range[end] * x_scale, 0)) 
    # closepath()
    
    # Fill
    setopacity(0.15)
    fillpreserve() 
    
    # Outline
    setopacity(1.0)
    setline(6)
    strokepath()   
    
    # -----------------------------------------------------
    # Layer 2: The Ordinate/Sigmoid (Green)
    # -----------------------------------------------------
    sethue(C_GREEN)
    setline(6)
    setlinecap("round")
    
    # Draw S-Curve shifted UP by vertical_offset
    move(Point(-3.0 * x_scale, -ordinate(-3.0) * y_scale - 0vertical_offset))
    for t in range(-3.0, 3.0, length=100)
        # We subtract vertical_offset to move it higher in the image
        line(Point(t * x_scale, -ordinate(t) * y_scale - 0vertical_offset))
    end
    strokepath()
    
    # -----------------------------------------------------
    # Layer 3: The Quadrature Nodes (Red)
    # -----------------------------------------------------
    sethue(C_RED)
    
    node_indices = -2:2
    
    for k in node_indices
        t_k = Float64(k) * 0.8
        
        # Calculate Position (Apply same shift)
        px = t_k * x_scale
        py = -ordinate(t_k) * y_scale - 0vertical_offset
        
        # Draw dot (Size reduced to 10, no rim)
        circle(Point(px, py), 8, :fill)
    end

    finish()
    println("Saved $filename")
end

draw_logo_centered("docs/FastTanhSinh_Logo.svg")