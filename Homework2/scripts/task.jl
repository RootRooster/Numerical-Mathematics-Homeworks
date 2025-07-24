using Homework2

CONTROL_POLYGON = [
  (0.0, 0.0),
  (1.0, 1.0),
  (2.0, 3.0),
  (1.0, 4.0),
  (0.0, 4.0),
  (-1.0, 3.0),
  (0.0, 1.0),
  (1.0, 0.0)
]

function main()
  println("\n--------------------")
  println("Calculating the area of the Bezier curve...")
  println("Control polygon:")
  for point in CONTROL_POLYGON
    println(" ", point)
  end

  area = calculate_bezier_loop_area(CONTROL_POLYGON)
  println("\n--------------------")
  println("Calculated area: ", area)
  println("--------------------")
end

main()

