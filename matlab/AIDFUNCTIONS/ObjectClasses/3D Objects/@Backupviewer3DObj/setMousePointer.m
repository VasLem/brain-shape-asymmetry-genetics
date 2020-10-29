function setMousePointer(obj)
  switch obj.Status
      case 'Ready'
        switch obj.ActiveKey
            case 'space'
                 switch obj.SelectionMode
                     case 'landmark'
                         setptr(obj.Figure,'fullcrosshair');
                     case {'area' 'fill' 'brush' 'full'}
                         setptr(obj.Figure,'arrow');
                     otherwise
                 end
            otherwise
                switch obj.Mode
                    case 'camera'
                        switch obj.Action
                            case 'none'
                                 setptr(obj.Figure,'hand');
                            case 'rotate camera'
                                 setptr(obj.Figure,'rotate');
                            case 'pan camera'
                                 setptr(obj.Figure,'closedhand');
                            case 'zoom camera'
                                setptr(obj.Figure,'glass');
                            case 'zoom in'
                                setptr(obj.Figure,'glassplus');
                            case 'zoom out'
                                setptr(obj.Figure,'glassminus');
                            otherwise
                                setptr(obj.Figure,'hand');
                        end
                    case 'light'
                        switch obj.Action
                            case 'none'
                                setptr(obj.Figure,'hand1');
                            case 'rotate light'
                                setptr(obj.Figure,'closedhand');
                            otherwise
                                setptr(obj.Figure,'hand1');
                        end                           
                    otherwise
                        setptr(obj.Figure,'arrow');
                end
        end
      case 'Busy'
          setptr(obj.Figure,'watch');
      otherwise
  end
end