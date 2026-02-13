from torch.optim.lr_scheduler import ReduceLROnPlateau

def train_one_epoch(model, loader, optimizer, criterion):
    model.train()
    total_loss = 0
    for data in loader:
        optimizer.zero_grad()
        output = model(data, data.target)
        loss = criterion(output, data.y)
        loss.backward()
        optimizer.step()
        total_loss += loss.item()
    return total_loss / len(loader)

# Optimization Strategy:
# - Adam Optimizer for stable convergence
# - ReduceLROnPlateau to fine-tune weights as loss plateaus
# - MSELoss as the standard for affinity regression
