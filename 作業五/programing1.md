**Ⅰ. Forward Euler Method**
```python
def forward_euler(h):
    t = np.arange(0, 10 + h, h)
    y = np.zeros(len(t))
    y[0] = 1
    for n in range(len(t) - 1):
        y[n+1] = y[n] * (1 - 5*h)
    return t, y
```
```python
plt.subplot(1,2,1)
for h in hs:
    t, y = forward_euler(h)
    plt.plot(t, y, 'o-', label=f'h={h}')
plt.plot(t_exact, y_forward_exact, 'k--', label='Exact')
plt.title('Forward Euler: y\' = -5y')
plt.xlabel('t'); plt.ylabel('y(t)')
plt.legend(); plt.grid(True)
```

**Forward Euler Result Analysis**

<table>
<tr>
<td width="45%">
<img src="圖一.png" width="100%">
</td>
<td width="55%">

### Explanation

- When  $$y' = -5y$$ , the analytical solution decays rapidly toward zero.
- The Forward Euler update rule is $$\[y_{n+1} = (1 - 5h) y_n\].$$
- If $$h < 0.4$$, the numerical solution remains stable.  
- When $$h = 0.41$$, we have  $$1 - 5h = -1.05$$.  
  The solution starts to oscillate and diverge, showing that the Forward Euler method becomes unstable for stiff problems.

</td>
</tr>
</table>


---

**Ⅱ. Backward Euler Method**

```python
def backward_euler(h):
    t = np.arange(0, 10 + h, h)
    y = np.zeros(len(t))
    y[0] = 1
    for n in range(len(t) - 1):
        y[n+1] = y[n] / (1 - 5*h)
    return t, y
```
```python
plt.subplot(1,2,2)
for h in hs:
    t, y = backward_euler(h)
    plt.plot(t, y, 'o-', label=f'h={h}')
plt.plot(t_exact, y_backward_exact, 'k--', label='Exact')
plt.title('Backward Euler: y\' = 5y')
plt.xlabel('t'); plt.ylabel('y(t)')
plt.legend(); plt.grid(True)
```
**Backward Euler Result Analysis**

<table>
<tr>
<td width="45%">
<img src="圖二.png" width="100%">
</td>
<td width="55%">

**Explanation**

- For $$y' = 5y$$, the analytical solution grows exponentially.
- The Backward Euler update rule is  $$\[y_{n+1} = \frac{y_n}{1 - 5h}\].$$
- Even for large step sizes \( h \), the numerical solution remains stable, although the error may increase.
- This demonstrates that the Backward Euler method is A-stable, making it suitable for stiff differential equations.

</td>
</tr>
</table>
