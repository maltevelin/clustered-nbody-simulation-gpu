--
-- ==
--
-- nobench input @ data/nbody-acc-t0.in
-- output @ data/nbody-acc-t0.out
--
-- nobench input @ data/nbody-acc-t10.in
-- output @ data/nbody-acc-t10.out
--
-- notest input @ data/1000-bodies.in
-- notest input @ data/10000-bodies.in
-- notest input @ data/100000-bodies.in
-- notest input @ data/200000-bodies.in
-- notest input @ data/1000000-bodies.in
-- notest input @ data/2000000-bodies.in
-- notest input @ data/4000000-bodies.in

import "lib/github.com/athas/matte/colour"
import "octree"

type acceleration = vec3.vector

let accel (epsilon: f32) (x:body) (y: {position: position, mass: mass})
  : velocity =
  let r = vec3.(y.position - x.position)
  let rsqr = vec3.dot r r + epsilon * epsilon
  let invr = 1.0f32 / f32.sqrt rsqr
  let invr3 = invr * invr * invr
  let s = y.mass * invr3
  in vec3.scale s r


let move [n] (epsilon: f32) (theta: f32)
             (t: [n]octnode, side_max: f32, root_delta: i32) (x: body)
             : acceleration =
  let (acceleration, _, _, _) =
  loop (acc, cur, prev, to_check) = ({x=0f32, y=0f32, z=0f32}, 0, -1, 1)
  while to_check > 0 do
    let from_child = prev != -1 && unsafe t[prev].parent == cur
    let current = unsafe t[cur]
    in
    if from_child then -- prev is my child => just move to next node
      --- index 'i' of child will be found before 'i' goes out of bounds.
      let j = loop i = 0 while (unsafe (current.children)[i]) != prev do
                         (i + 1)
      let next_cur =
        if j == 7 || (unsafe (current.children)[j + 1]) == -1
          then current.parent
        else unsafe (current.children)[j + 1] --- j must be less than 7
      in (acc, next_cur, cur, to_check)
    else -- prev is my parent
      let s = side_max / 2.0f32**(f32.i32 current.tree_level)
      let position = current.body.position
      let total_mass = current.body.mass
      let q = vec3.scale (1.0f32 / total_mass) position
      let p = x.position
      let inner =
        (p.x - q.x)**2.0f32 + (p.y - q.y)**2.0f32 + (p.z - q.z)**2.0f32
      let local_theta =  s / f32.sqrt inner
      let overlapping_subtree = (root_delta / 3) + current.tree_level > 11
      in
      -- node is far enough, so we update
      if current.is_leaf || overlapping_subtree || local_theta < theta then
        let velocity = accel epsilon x {position = q, mass = total_mass}
        in (acc vec3.+ velocity, current.parent, cur, to_check - 1)
      else -- node is too close
        let children = current.children
        let num_children = map (i32.bool <-< (>= 0)) children
                           |> reduce_comm (+) 0
        in (acc, children[0], cur, to_check + num_children - 1)
  in acceleration

let calc_accels [n] (epsilon: f32) (theta: f32) (bodies: [n]body)
                    : ([n]acceleration, [n]body) =
  let (accelerator, side_max, root_delta, sorted_bodies) =
    mk_accelerator bodies
  in (map (move epsilon theta (accelerator, side_max, root_delta)
      ) sorted_bodies, sorted_bodies)


let advance_body (time_step: f32) (body: body) (acc: acceleration): body =
  let acceleration = vec3.scale body.mass acc
  let position = vec3.(body.position + scale time_step body.velocity)
  let velocity = vec3.(body.velocity + scale time_step acceleration)
  in {position, mass=body.mass, velocity}

let advance_bodies [n] (epsilon: f32) (time_step: f32) (theta: f32)
                       (bodies: [n]body): [n]body =
  let (accels, bodies) = calc_accels epsilon theta bodies
  in map2 (advance_body time_step) bodies accels

let advance_bodies_steps [n] (n_steps: i32) (epsilon: f32) (time_step: f32)
                             (theta: f32) (bodies: [n]body): [n]body =
  loop bodies for _i < n_steps do
    advance_bodies epsilon time_step theta bodies

let wrap_body (posx: f32, posy: f32, posz: f32)
              (mass: f32)
              (velx: f32, vely: f32, velz: f32): body =
  {position={x=posx, y=posy, z=posz},
   mass,
   velocity={x=velx, y=vely, z=velz}}

let unwrap_body (b: body): ((f32, f32, f32), f32, (f32, f32, f32)) =
  ((b.position.x, b.position.y, b.position.z),
   b.mass,
   (b.velocity.x, b.velocity.y, b.velocity.z))

entry main [n]
        (n_steps: i32)
        (epsilon: f32)
        (time_step: f32)
        (theta: f32)
        (xps: [n]f32)
        (yps: [n]f32)
        (zps: [n]f32)
        (ms: [n]f32)
        (xvs: [n]f32)
        (yvs: [n]f32)
        (zvs: [n]f32): ([n]f32, [n]f32, [n]f32, [n]f32, [n]f32, [n]f32, [n]f32) =
  let bodies  = map3 wrap_body (zip3 xps yps zps) ms (zip3 xvs yvs zvs)
  let bodies' = advance_bodies_steps n_steps epsilon time_step theta bodies
  let (final_pos, ms', final_vel) = map unwrap_body (bodies') |> unzip3
  let (xps', yps', zps') = unzip3 final_pos
  let (xvs', yvs', zvs') = unzip3 final_vel
  in (xps', yps', zps', ms', xvs', yvs', zvs')

let rotatePointByMatrix (rotation: [3][3]f32) ({x,y,z}: position): position =
  {x= x*rotation[0,0] + y*rotation[1,0] + z*rotation[2,0],
   y= x*rotation[0,1] + y*rotation[1,1] + z*rotation[2,1],
   z= x*rotation[0,2] + y*rotation[1,2] + z*rotation[2,2]}

let rotatePointsByMatrix [n] (rotation: [3][3]f32)(ps: [n]position): [n]position =
  map (rotatePointByMatrix rotation) ps

let rotateXMatrix (angle: f32): [3][3]f32 =
  [[1f32,        0f32,         0f32],
   [0f32, f32.cos angle, -f32.sin angle],
   [0f32, f32.sin angle,  f32.cos angle]]

let rotateYMatrix (angle: f32): [3][3]f32 =
  [[f32.cos angle,  0f32, f32.sin angle],
   [0f32,           1f32, 0f32],
   [-f32.sin angle, 0f32, f32.cos angle]]

let matmult [n][m][p] (x: [n][m]f32) (y: [m][p]f32): [n][p]f32 =
  map (\xr ->
         map (\yc -> f32.sum (map2 (*) xr yc))
             (transpose y))
      x

let rotationMatrix (x_rotation: f32) (y_rotation: f32): [3][3]f32 =
  matmult (rotateXMatrix x_rotation) (rotateYMatrix y_rotation)

let inverseRotationMatrix (x_rotation: f32) (y_rotation: f32): [3][3]f32 =
  matmult (rotateYMatrix y_rotation) (rotateXMatrix x_rotation)

entry inverseRotatePoint (x: f32) (y: f32) (z: f32) (x_rotation: f32) (y_rotation: f32): (f32,f32,f32) =
  let {x,z,y} = rotatePointByMatrix (inverseRotationMatrix x_rotation y_rotation) {x,y,z}
  in (x,y,z)

let rotatePoints [n] (ps: [n]position) (x_rotation: f32) (y_rotation: f32): [n]position =
  rotatePointsByMatrix (rotationMatrix x_rotation y_rotation) ps

let renderPoint(h: i32) (w: i32) (x_ul: f32) (y_ul: f32) (x_br: f32) (y_br: f32) (max_mass: f32)
               ({x,y,z=_}:position) (m: f32): (i32, i32) =
  -- Draw nothing if the point is outside the viewport.
  if x < x_ul || x > x_br || y < y_ul || y > y_br then (-1, 0)
  else
  -- Normalise x,y to positions in interval (0,1) within the viewport.
  let x' = (x-x_ul) / (x_br-x_ul)
  let y' = (y-y_ul) / (y_br-y_ul)
  -- Convert x',y' to screen coordinate space.
  let x'' = t32(x' * r32(w))
  let y'' = t32(y' * r32(h))
  let intensity = if m >= max_mass
                  then 255
                  else 128 + t32((m / max_mass) * 128f32)
  let colour = intensity * 0x10000 +
               intensity * 0x100 +
               0xFF
  in (y''*w + x'', colour)

entry render [n]
            (h: i32) (w: i32) (x_ul: f32) (y_ul: f32) (x_br: f32) (y_br: f32)
            (xps: [n]f32) (yps: [n]f32) (zps: [n]f32) (ms: [n]f32)
            (x_rotation: f32) (y_rotation: f32)
            (max_mass: f32) (invert: bool): [h][w]i32 =
  let background = if invert then argb.white else argb.black
  let (is, vs) = unzip(map2 (renderPoint h w x_ul y_ul x_br y_br max_mass)
                            (rotatePoints (map3 (\x y z -> {x,y,z}) xps yps zps)
                                          x_rotation y_rotation) ms)
  let vs' = map (\x -> if invert then !x else x) vs
  in unflatten h w (scatter (replicate (h*w) background) is vs')

entry mouse_mass_active (xps: *[]f32) (yps: *[]f32) (zps: *[]f32) (ms: *[]f32) (x: f32) (y: f32) (z: f32) =
  (xps with [0] = x,
   yps with [0] = y,
   zps with [0] = z,
   ms with [0] = 10000)

entry mouse_mass_inactive (xps: *[]f32) (yps: *[]f32) (zps: *[]f32) (ms: *[]f32) =
  (xps with [0] = 0,
   yps with [0] = 0,
   zps with [0] = 0,
   ms with [0] = 0.0001)
