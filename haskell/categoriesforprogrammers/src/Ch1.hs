module Ch1 where

import Data.Void as V
import Data.Serialize as S
import Data.Map as M
import Data.ByteString as B
import System.Directory as D


f :: Bool -> Bool
f = undefined -- f maps any bool to _|_, it is also element (_|_) of Bool -> Bool

fact n = product [1..n]

absurd :: V.Void -> a
absurd = undefined -- cannot be called, because there are no values of type Void

f44 :: () -> Integer
f44 () = 44 -- functions from singleton set are tantamount to selecting elements from sets

fInt :: Integer -> ()
fInt _ = () -- example function to singleton set

unit :: a -> ()
unit _ = () -- example of parametric polymorphism

memoFile = "/Users/kbrown/scripts/haskell/categoriesforprogrammers/src/Ch1.txt"

memoizeFact :: Integer ->  IO ()
memoizeFact a = do 
                memoexists <- D.doesPathExist memoFile
                if memoexists 
                    then 
                        do 
                            memostr <- B.readFile memoFile
                            let memo = case S.decode memostr of 
                                            Left _ -> M.empty
                                            Right m ->  m 

                            let known =  M.member a memo
                            if known
                                then print $ memo ! a
                                else writenew a (fact a) memo
                    else 
                        writenew a (fact a) M.empty

    where writenew k v mp = do  let newmemo = M.insert k v mp 
                                B.writeFile memoFile (S.encode newmemo) 
                                print v

true :: Bool -> Bool
true _ = True
false :: Bool -> Bool
false _ = False
idBool :: Bool -> Bool
idBool = id
notBool :: Bool -> Bool
notBool = not