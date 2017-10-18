module Misc where
import Linear 			as L
import Linear.Matrix 	as M
-- import etc.

if' :: Bool -> a -> a -> a
if' True  x _ = x
if' False _ y = y

cell2param :: M.M33 Double -> L.V3 Double
cell2param = undefined

symbols2electrons :: [String] -> Int
symbols2electrons = undefined