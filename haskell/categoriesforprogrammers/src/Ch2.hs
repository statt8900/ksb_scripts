module Ch2 where
import Data.Map as M
map' = M.map.fromList
data Category = Cat {c_ob   :: Integer
					,c_hom  :: Integer
					,c_src  :: M.Map Integer Integer
					,c_tgt  :: M.Map Integer Integer
					,c_comp :: M.Map (Integer,Integer) Integer} 
-- specifies a FINITE category up to isomorphism


data PreOrder = PrO {pro_ob :: Integer
                    ,pro_rel:: (Integer,Integer)}
data ParOrder = ParO {par_ob :: Integer
                    ,par_rel:: (Integer,Integer)}
data LinOrder = LinO {lin_ob :: Integer
                    ,lin_rel:: (Integer,Integer)}


type Funct = Category -> Category

zero = Cat(0,0,M.empty,M.empty,M.empty)
one = Cat(1,1,map' [(1,1)],map' [(1,1)])

